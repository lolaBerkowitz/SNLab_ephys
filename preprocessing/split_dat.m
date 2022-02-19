function split_dat(file_name, folder_order,varargin)
% split_dat for recordings of multiple animals from Intan RHX USB
% acquisition system. Assumes intan one data per data type format. Will
% split amplifer and aux by input ports by default.


% input:
%   file_location: name of file location with data
%   folder_order: cell array with folder location (i.e. {'HPC01','HPC02'...}
%            will give data_path\HPC01, data_path\HPC02, ...). Order must
%            match port_order (default is ABCD)
%
% variable arguments:
%   split_folder: string indicating folder where session files for headstage
%           recordings are kept. Should be within the main Data folder

%           where folder locations are destined i.e. data_path\to_split.
%   digitalin_order: numeric vector that signifies the digitalin channel
%           that corresponds to the port_order. Default is [2,3,4,5] so
%           events for port A would be on channel 1, port B channel 2, etc.
%
%   port_order: cell array containing order of ports consistent with folder order, so {'A','B',...}
%            would result with folder HPC01 having channels saved to dat from Port A, HPC02
%            having channels saved to dat from Port B, etc.

% TO-DO:
%  - Make dynamic input so multiple ports can be assiged to the same animal
%    for drives over 64 channels.
% - Add availablity to scan data_path\split_folder and automatically process folders that
%    have not yet been processed


% output:
%   saves the dat file split by occupied ports into folders.

p = inputParser;
addParameter(p,'split_folder','to_split',@isstring)
addParameter(p,'data_path','D:\app_ps1\data',@isstring)
addParameter(p,'digitalin_order',[2,3,4,5]',@isnumeric) %snlab rig wiring
addParameter(p,'port_order',{'A','B','C','D'},@isarray)
parse(p,varargin{:});

split_folder = p.Results.split_folder;
data_path = p.Results.data_path;
digitalin_order = p.Results.digitalin_order;
port_order = p.Results.port_order;

dat_folder = [data_path,filesep,split_folder,filesep,file_name];

% Load info.rhd for port information
[amplifier_channels, ~, aux_input_channels, ~,...
    ~, ~, frequency_parameters,~ ] = ...
    read_Intan_RHD2000_file_snlab(dat_folder);

% splits dat according to folder_order and saves to dat_folder
% process_aux(dat_folder,aux_input_channels,frequency_parameters,folder_order);
% process_amp(dat_folder,amplifier_channels,frequency_parameters,folder_order)

% load digitalin channels for session start end times (first and last event per channel)
digitalIn = process_digitalin(dat_folder,'digitalin.dat',frequency_parameters.board_dig_in_sample_rate);

%% Loop through folders.
% Order of folders should match port_order and digitalin_order
for i = find(~cellfun(@isempty,folder_order))
    
    % make basepath for folders
    basepath = [data_path,filesep,folder_order{i},filesep,[folder_order{i},'_',file_name]];
    
    if ~exist(basepath,'dir')
        mkdir(basepath)
    end
    % move amplifier and aux files to basepath and rename to standard intan
    % format
    if isfile([dat_folder,filesep,port_order{i}(end),'_amplifier.dat'])
        movefile([dat_folder,filesep,port_order{i}(end),'_amplifier.dat'],basepath);
        system(['rename ',basepath,filesep,port_order{i}(end),'_amplifier.dat',' ','amplifier.dat']);
    end
    if isfile([dat_folder,filesep,port_order{i}(end),'_aux.dat'])
        movefile([dat_folder,filesep,port_order{i}(end),'_aux.dat'],basepath);
        system(['rename ',basepath,filesep,port_order{i}(end),'_aux.dat',' ','auxillary.dat']);
    end
    % if video file with subid is present, move that to basepath
    if isfile([dat_folder,filesep,'*_',folder_order{i},'.avi'])
        movefile([dat_folder,filesep,'*_',folder_order{i},'.avi'],basepath)
    end
    
    % make copy of rhd, setting, and time to basepath
    copyfile([dat_folder,filesep,'time.dat'],basepath);
    copyfile([dat_folder,filesep,'settings.xml'],basepath);
    copyfile([dat_folder,filesep,'info.rhd'],basepath);
    
    % create digitalIn event structure and save to basepath
    parse_digitalIn(digitalIn,digitalin_order(i),basepath)
    
end


end

% function main(basepath,amp,aux,time)


function process_aux(dat_path,aux_input_channels,frequency_parameters,ports)
% processes intan auxillary file and splits into separate files indicated
% by ports.
% input:
% - dat_path: path where auxiliary.dat is found
% - aux_input_channels: table containing metadata from intan RHD file
%   (output of read_Intan_RHD2000_file_snlab)
% - frequency_parameters: table containing metadata from intan RHD file
%   (output of read_Intan_RHD2000_file_snlab)
% - ports: cell array containing
% output:
% - saves port_aux.dat to dat_path

% Check to see if files have been created 
% find ports to write
write_port = unique([{aux_input_channels.port_name}]);
write_port = write_port(find(~cellfun(@isempty,ports))); % write ports based on inputs

% check if files have been created and if so update list
for port = write_port
    if isfile([dat_path,filesep,port{1}(end),'_aux.dat'])
        disp([port{1}(end),'_aux.dat ','already created'])
        % remove from list
        write_port(find(ismember(write_port,port))) = [];
    end
end

% if all files have been written, exit the function
if isempty(write_port)
    return
end

% Load the file
n_channels = size(aux_input_channels,2);
contFile = fullfile(dat_path,'auxiliary.dat');
file = dir(contFile);
samples = file.bytes/(n_channels * 2); %int16 = 2 bytes
aux = memmapfile(contFile,'Format',{'uint16' [n_channels samples] 'mapped'});
x = aux.Data.mapped;

% create batches
batch = ceil(linspace(0,samples,ceil(samples/frequency_parameters.aux_input_sample_rate/4)));

% loop through ports
for port = write_port
    
    % skip if file is already created
    if isfile([dat_path,filesep,port{1}(end),'_aux.dat'])
        disp([port{1}(end),'_aux.dat ','already created'])
        continue
    end
    
    % initiate file
    aux_file = fopen(fullfile(dat_path,[port{1}(end),'_aux' ,'.dat']),'w');
    idx = contains({aux_input_channels.port_name},port{1});
    
    % loop though batches save port channels to file
    for i = 1:length(batch)-1
        disp(['batch ',num2str(batch(i)+1),' to ',num2str(batch(i+1)),...
            '   ',num2str(i),' of ',num2str(length(batch)-1)])
        
        % pull out batch of data
        datr = x(idx,batch(i)+1:batch(i+1)) * 0.0000374; % 0.0000374 = conver to volts
        
        % gather data and convert to int16 & transpose
        datcpu = gather_try(int16(datr));
        
        % write to disk
        fwrite(aux_file, datcpu, 'int16');
    end
    fclose(aux_file);
end
clear x aux 
end

function process_amp(dat_path,amplifier_channels,frequency_parameters,ports)
% processes intan amplifier file and splits into separate files indicated
% by ports.
% input:
% - dat_path: path where amplifier.dat is found
% - amplifier_channels: table containing metadata from intan RHD file
%   (output of read_Intan_RHD2000_file_snlab)
% - frequency_parameters: table containing metadata from intan RHD file
%   (output of read_Intan_RHD2000_file_snlab)
% - ports: cell array containing
% output:
% - saves port_amplifier.dat to dat_path

% Check to see if files have been created 
% loop through ports
write_port = unique([{amplifier_channels.port_name}]);
write_port = write_port(find(~cellfun(@isempty,ports))); % write ports based on inputs

% loop through ports
for port = write_port
    if isfile([dat_path,filesep,port{1}(end),'_amplifier.dat'])
        disp([port{1}(end),'_amplifier.dat ','already created'])
        % remove from list
        write_port(find(ismember(write_port,port))) = [];
    end
end

% if all files have been written, exit the function
if isempty(write_port)
    return
end

n_channels = size(amplifier_channels,2);
contFile = fullfile(dat_path,'amplifier.dat');
file = dir(contFile);
samples = file.bytes/(n_channels * 2); %int16 = 2 bytes
amp = memmapfile(contFile,'Format',{'int16' [n_channels samples] 'mapped'});
% amp.Data = amp.Data.Data.mapped * 0.195; % covert to microvolts

x = amp.Data.mapped;
% aux.Data = aux.Data * 0.0000374; % covert to volts

% create batches
batch = ceil(linspace(0,samples,ceil(samples/frequency_parameters.amplifier_sample_rate/4)));

% loop through ports
for port = write_port
    
    amp_file = fopen(fullfile(dat_path,[port{1}(end),'_amplifier' ,'.dat']),'w');
    idx = contains({amplifier_channels.port_name},port{1});
    
    % loop though batches
    for i = 1:length(batch)-1
        disp(['batch ',num2str(batch(i)+1),' to ',num2str(batch(i+1)),...
            '   ',num2str(i),' of ',num2str(length(batch)-1)])
        
        % pull out batch of data
        datr = x(idx,batch(i)+1:batch(i+1))* 0.195; % 0.195 = convert to uV
        
        % gather data and convert to int16 & transpose
        datcpu = gather_try(int16(datr));
        
        % write to disk
        fwrite(amp_file, datcpu, 'int16');
    end
    
    fclose(amp_file);
end
clear x amp 
end

function digitalIn = process_digitalin(data_path,dat_name,fs)
% code adapted from getDigitalin.m in neurocode
% (https://github.com/ayalab1/neurocode/tree/master/preprocessing)

lag = 20; % This pertains to a period for known event (in this case stimulation period).
%           keeping for now but will need to update if we ever have events
%           for stimulation. Making large cause events are based on double-clicks
%           of varying length LB 2/22


contFile = fullfile(data_path,dat_name);
% file = dir(fullfile(data_path,dat_name));
% samples = file.bytes/2/16; % 16 is n_channels for intan digitial in
D.Data = memmapfile(contFile,'Format','uint16','writable',false);

digital_word2 = double(D.Data.Data);
Nchan = 16;
Nchan2 = 17;
for k = 1:Nchan
    tester(:,Nchan2-k) = (digital_word2 - 2^(Nchan-k))>=0;
    digital_word2 = digital_word2 - tester(:,Nchan2-k)*2^(Nchan-k);
    test = tester(:,Nchan2-k) == 1;
    test2 = diff(test);
    pulses{Nchan2-k} = find(test2 == 1);
    pulses2{Nchan2-k} = find(test2 == -1);
    data(k,:) = test;
end
digital_on = pulses;
digital_off = pulses2;


for ii = 1:size(digital_on,2)
    if ~isempty(digital_on{ii})
        % take timestamp in seconds
        digitalIn.timestampsOn{ii} = digital_on{ii}/fs;
        digitalIn.timestampsOff{ii} = digital_off{ii}/fs;
        
        % intervals
        d = zeros(2,max([size(digitalIn.timestampsOn{ii},1) size(digitalIn.timestampsOff{ii},1)]));
        d(1,1:size(digitalIn.timestampsOn{ii},1)) = digitalIn.timestampsOn{ii};
        d(2,1:size(digitalIn.timestampsOff{ii},1)) = digitalIn.timestampsOff{ii};
        if d(1,1) > d(2,1)
            d = flip(d,1);
        end
        if d(2,end) == 0; d(2,end) = nan; end
        digitalIn.ints{ii} = d;
        digitalIn.dur{ii} = digitalIn.ints{ii}(2,:) - digitalIn.ints{ii}(1,:); % durantion
        
        clear intsPeriods
        intsPeriods(1,1) = d(1,1); % find stimulation intervals
        intPeaks =find(diff(d(1,:))>lag);
        for jj = 1:length(intPeaks)
            intsPeriods(jj,2) = d(2,intPeaks(jj));
            intsPeriods(jj+1,1) = d(1,intPeaks(jj)+1);
        end
        intsPeriods(end,2) = d(2,end);
        digitalIn.intsPeriods{ii} = intsPeriods;
    end
end
end

function parsed_digitalIn = parse_digitalIn(digitalIn,channel_index,basepath,varargin)
% saves events from digitalIn structure. By default, the video timestamp
% data is obtained from intan digitalin channel 0 on the RHD USB interface
% board
% input:
%  - digitalIn: structure produced by process_digitalIn containing events
%     from intan digitalIn.dat file
%  - channel_index: channels to be saved. Assum
p = inputParser;
addParameter(p,'video_idx',1,@isnumeric)
parse(p,varargin{:});

video_idx = p.Results.video_idx;

parsed_digitalIn.timestampsOn{1,video_idx} = digitalIn.timestampsOn{1, video_idx};
parsed_digitalIn.timestampsOff{1,video_idx} = digitalIn.timestampsOff{1, video_idx};
parsed_digitalIn.timestampsOn{1,2} = digitalIn.timestampsOn{1, channel_index};
parsed_digitalIn.timestampsOff{1,2} = digitalIn.timestampsOff{1, channel_index};
parsed_digitalIn.ints{1,2} = digitalIn.ints{1, channel_index};
parsed_digitalIn.dur{1,2} = digitalIn.dur{1, channel_index};
parsed_digitalIn.intsPeriods{1,2} = digitalIn.intsPeriods{1, channel_index};

save([basepath,filesep,'digitalIn.events.mat'],'parsed_digitalIn');

end

