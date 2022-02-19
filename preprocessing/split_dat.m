function split_dat(file_name, folder_order,varargin)
% split_dat for recordings of multiple animals from Intan RHX USB
% acquisition system. Assumes intan one data per data type format. Will  
% split amplifer and aux by available ports by default. 
% Optional: Uses assigned digital input channels to split the file at specific 
% time epochs. 

% input: 
%   file_location: name of file location with data
%   folder_order: cell array with folder location (i.e. {'HPC01','HPC02'...} 
%            will give data_path\HPC01, data_path\HPC02, ...). Order must
%            match port_order (default is ABCD) 
%
% variable arguments: 
%   split_folder: name of folder where session files for headstage
%           recordings are kept. Should be within the main Data folder 

%           where folder locations are destined i.e. data_path\to_split.
%   digitalin_order: numeric vector that signifies the digitalin channel 
%           that corresponds to the port_order. Default is [2,3,4,5] so 
%           events for port A would be on channel 1, port B channel 2, etc.            
%
%   port_order: order of ports consistent with folder order, so {'A','B',...}
%            would result with folder HPC01 having channels saved to dat from Port A, HPC02 
%            having channels saved to dat from Port B, etc.
% TO-DO: 
%  - Make dynamic input so multiple ports can be assiged to the same animal
%    for drives over 64 channels. 
% - Add availablity to scan drive and automatically process folders that
%    have not yet been processed

             
% output: 
%   saves the dat file split by occupied ports into folders. 

p = inputParser;
addParameter(p,'split_folder','to_split',@isstring) 
addParameter(p,'data_path','D:\app_ps1\data',@isstring)
addParameter(p,'digitalin_order',[2,3,4,5]',@isnumeric) 
addParameter(p,'port_order',{'A','B','C','D'},@isarray) 
parse(p,varargin{:});

split_folder = p.Results.split_folder;
data_path = p.Results.data_path;
digitalin_order = p.Results.digitalin_order;
port_order = p.Results.port_order;

% Load info.rhd for port information
dat_folder = [data_path,filesep,split_folder,filesep,file_name];

[amplifier_channels, ~, aux_input_channels, ~,...
    ~, ~, frequency_parameters,~ ] = ...
    read_Intan_RHD2000_file_snlab(dat_folder);
fs = frequency_parameters.amplifier_sample_rate;  
amp_channels = 1:size(amplifier_channels,2);
aux_channels = 1:size(aux_input_channels,2);


%% Load dat files
process_aux(dat_folder,aux_input_channels,frequency_parameters,folder_order);
process_amp(dat_folder,amplifier_channels,frequency_parameters,folder_order)

% load digitalin channels for session start end times (first and last event per channel)
digitalIn = process_digitalin(dat_folder,'digitalin.dat',frequency_parameters.board_dig_in_sample_rate);
time = load_time(dat_folder);

%% Loop through folders. 
% Order of folders should match port_order and digitalin_order 
for i = find(~cellfun(@isempty,folder_order))

    % find start and end of events for iteration 
    session_idx = [min(digitalIn.timestampsOn{1,digitalin_order(i)}*fs),max(digitalIn.timestampsOff{1,digitalin_order(i)}*fs)];
    basepath = [data_path,filesep,folder_order{i},filesep,[folder_order{i},'_',file_name]];

    switch port_order{i}

        case 'A'
            % find channels for port
            amp_index = amp_channels(contains({amplifier_channels.port_name},'Port A'));
            aux_index = aux_channels(contains({aux_input_channels.port_name},'Port A'));
            
        case 'B'
            amp_index = amp_channels(contains({amplifier_channels.port_name},'Port B'));
            aux_index = aux_channels(contains({aux_input_channels.port_name},'Port B'));
            
        case 'C'
            amp_index = amp_channels(contains({amplifier_channels.port_name},'Port C'));
            aux_index = aux_channels(contains({aux_input_channels.port_name},'Port C'));

        case 'D'
            amp_index = amp_channels(contains({amplifier_channels.port_name},'Port D'));
            aux_index = aux_channels(contains({aux_input_channels.port_name},'Port D'));
    end


    % index dat and write amplifier.dat and auxillary.dat to basepath
%     write_aux_dat(Aux,basepath,aux_index,session_idx)
% 
%     write_amp_dat(Amp,basepath,amp_index,session_idx)

    % create event file
    parse_digitalIn(digitalIn,session_idx,digitalin_order(i),basepath)
            
end

% Remove signal outside of start and end event epochs 

% Save dat to animal folder 


end

% function main(basepath,amp,aux,time)


function process_aux(dat_path,aux_input_channels,frequency_parameters,ports)

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
write_port = unique([{aux_input_channels.port_name}]);
write_port = write_port(find(~cellfun(@isempty,ports))); % write ports based on inputs

for port = write_port 
    
    % skip if file is already created
    if isfile([dat_path,filesep,port{1}(end),'_aux.dat'])
        disp([port{1}(end),'_aux.dat_','already created'])
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
        datr = x(idx,batch(i)+1:batch(i+1));

        % gather data and convert to int16 & transpose
        datcpu = gather_try(int16(datr));

        % write to disk
        fwrite(aux_file, datcpu, 'int16');
    end
    fclose(aux_file);
end
end

function process_amp(dat_path,amplifier_channels,frequency_parameters,ports)

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
write_port = unique([{amplifier_channels.port_name}]);
write_port = write_port(find(~cellfun(@isempty,ports))); % write ports based on inputs

% loop through ports
for port = write_port
    if isfile([dat_path,filesep,port{1}(end),'_amplifier.dat'])
        continue
    end
    
    amp_file = fopen(fullfile(dat_path,[port{1}(end),'_amplifier' ,'.dat']),'w');
    idx = contains({amplifier_channels.port_name},port{1});
    
    % loop though batches
    for i = 1:length(batch)-1
        disp(['batch ',num2str(batch(i)+1),' to ',num2str(batch(i+1)),...
            '   ',num2str(i),' of ',num2str(length(batch)-1)])

        % pull out batch of data
        datr = x(idx,batch(i)+1:batch(i+1));

        % gather data and convert to int16 & transpose
        datcpu = gather_try(int16(datr));

        % write to disk
        fwrite(amp_file, datcpu, 'int16');
    end
    fclose(aux_file);
end
end

function time = load_time(dat_path)
contFile = fullfile(dat_path,'time.dat');
fileinfo = dir(contFile);
num_samples = fileinfo.bytes/4; % int32 = 4 bytes
time.Data = memmapfile(contFile,'Format',{'int32' [1 num_samples] 'mapped'});
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

function digitalIn_new = parse_digitalIn(digitalIn,session_idx,channel_index)

digitalIn_new.timestampsOn{1,1} = digitalIn.timestampsOn{1, 1};
digitalIn_new.timestampsOff{1,1} = digitalIn.timestampsOff{1, 1};
digitalIn_new.timestampsOn{1,2} = digitalIn.timestampsOn{1, channel_index};
digitalIn_new.timestampsOff{1,2} = digitalIn.timestampsOff{1, channel_index};

end

function save_events(basepath, parsed_digitalIn)
% save events to basepath
save([basepath,filesep,'digitalIn.events.mat'],'parsed_digitalIn');
end
