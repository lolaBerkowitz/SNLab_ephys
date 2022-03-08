function split_dat(data_path,save_path, subject_order,varargin)
% split_dat for recordings of multiple animals from Intan RHX USB
% acquisition system. Assumes intan one data per data type format. Will
% split amplifer and aux by input ports by default.


% input:
%   data_path: path to file location with intan data to be split
%   save_path: cell array of paths to general directory where subject data folders are
%              located (one path per port)
%   subject_order: cell array with subject name (i.e. {'HPC01','HPC02'...}
%            will give data_path\HPC01, data_path\HPC02, ...). Order must
%            match order: Port A, Port B, Port C, Port D i.e. HPC01 was recording on 
%            Port A, HPC02 on Port B etc.
%
% variable arguments:
%   digitalin_order: numeric vector that signifies the digitalin channel
%           that corresponds to the port_order. Default is [2,3,4,5] so
%           events for port A would be on channel 1, port B channel 2, etc.
%

% TO-DO:
%  - Make dynamic input so multiple ports can be assiged to the same animal
%    for drives over 64 channels.
% - Add availablity to scan data_path\split_folder and automatically process folders that
%    have not yet been processed



% output:
%   saves the dat file split by occupied ports into folders.
%
% L Berkowitz 03/2022
p = inputParser;
addParameter(p,'digitalin_order',[2,3,4,5]',@isnumeric) %snlab rig wiring
parse(p,varargin{:});

digitalin_order = p.Results.digitalin_order;

% locate folder 
file_name = basenameFromBasepath(data_path);

% Load info.rhd for port information
[amplifier_channels, ~, aux_input_channels, ~,...
    ~, ~, frequency_parameters,~ ] = ...
    read_Intan_RHD2000_file_snlab(data_path);

% make the basepath if its not already made
for i = find(~cellfun(@isempty,subject_order))
    basepath{i} = fullfile(save_path{i},filesep,subject_order{i},[subject_order{i},'_',file_name]);
    
    if ~isfolder(basepath{i})
        mkdir(basepath{i})
    end
end

% splits dat according to subject_order and saves to data_path
process_aux(data_path,aux_input_channels,subject_order,basepath);
process_amp(data_path,amplifier_channels,frequency_parameters,subject_order,basepath);

% load digitalin channels for session start end times (first and last event per channel)
digitalIn = process_digitalin(data_path,'digitalin.dat',frequency_parameters.board_dig_in_sample_rate);

%% Loop through folders.
% Order of folders should match port_order and digitalin_order
for i = find(~cellfun(@isempty,subject_order))
    
    % if video file with subid is present, move that to basepath
    if isfile([data_path,filesep,'*_',subject_order{i},'.avi'])
        movefile([data_path,filesep,'*_',subject_order{i},'.avi'],basepath{i})
    end
    
    % make copy of rhd, setting, and time to basepath
    copyfile([data_path,filesep,'time.dat'],basepath{i});
    copyfile([data_path,filesep,'settings.xml'],basepath{i});
    copyfile([data_path,filesep,'info.rhd'],basepath{i});
    
    % create digitalIn event structure and save to basepath
    parse_digitalIn(digitalIn,digitalin_order(i),basepath{i})
end


end

% function main(basepath,amp,aux,time)


function process_aux(dat_path,aux_input_channels,folders,basepath)
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
% LB 03/22

% Check to see if files have been created 
% find ports to write
write_port = unique({aux_input_channels.port_name});
basepath = basepath(~cellfun(@isempty,folders));
write_port = write_port(~cellfun(@isempty,folders)); % write ports based on inputs

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

% loop through ports
tic
parfor port = 1:length(write_port)
    process_aux_(aux,write_port{port},aux_input_channels,basepath{port})
end
toc
clear aux 
end

function process_aux_(aux,port,aux_input_channels,basepath)
% LB 03/22
    % skip if file is already created
    if isfile([basepath,filesep,'auxiliary.dat'])
        disp([basepath,'auxiliary.dat ','already created'])
        return
    end
    
    % initiate file
    aux_file = fopen([basepath,filesep,'auxiliary.dat'],'w');
    idx = contains({aux_input_channels.port_name},port);
   
    % write to disk
    fwrite(aux_file, aux.Data.mapped(idx,:), 'uint16');
    fclose(aux_file);
end

function process_amp(dat_path,amplifier_channels,frequency_parameters,subject_order,basepath)
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
% LB 03/22
% Check to see if files have been created 
% loop through ports
write_port = unique({amplifier_channels.port_name});
basepath = basepath(~cellfun(@isempty,subject_order));
write_port = write_port(~cellfun(@isempty,subject_order)); % write ports based on inputs

% loop through ports
remove_port_idx = [];
for port = 1:length(write_port)
    if isfile([basepath{port},filesep,'amplifier.dat'])
        disp([basepath{port},' amplifier.dat ','already created'])
        % remove from list
        remove_port_idx = [remove_port_idx;find(ismember(write_port,write_port{port}))];
    end
end

if ~isempty(remove_port_idx)
    basepath(remove_port_idx) = [];
    write_port(remove_port_idx) = [];
end

% if all files have been written, exit the function
if isempty(write_port)
    return
end

n_channels = size(amplifier_channels,2);
contFile = fullfile(dat_path,'amplifier.dat');
file = dir(contFile);
samples = file.bytes/(n_channels * 2); %int16 = 2 bytes
amp = memmapfile(contFile,'Format',{'int16' [n_channels, samples] 'mapped'});

% create batches
batch = ceil(linspace(0,samples,ceil(samples/frequency_parameters.amplifier_sample_rate/4)));

% loop through ports
for port = 1:length(write_port)
    amp_file{port} = fopen(fullfile(basepath{port},'amplifier.dat'),'w');
    idx{port} = contains({amplifier_channels.port_name},write_port{port});
end
tic
% loop though batches
for i = 1:length(batch)-1
    disp(['batch ',num2str(batch(i)+1),' to ',num2str(batch(i+1)),...
        '   ',num2str(i),' of ',num2str(length(batch)-1)])
    % write to disk
    for port = 1:length(write_port)
        fwrite(amp_file{port},amp.Data.mapped(idx{port},batch(i)+1:batch(i+1)) * 0.195, 'int16');
    end
end

for port = 1:length(write_port)
    fclose(amp_file{port});
end
toc
clear amp 
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
%  - channel_index: channels to be saved. 
% LB 03/22
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

save([basepath,filesep,'digitalIn.events.mat'],'digitalIn');

end


