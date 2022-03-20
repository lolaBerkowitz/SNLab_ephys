function trim_dat(basepath,varargin)
% trim_dat removes parts of the dat file that does not contain neural data, 
% either from unplugs.

% Dependenies: one signal per channel dat files from intan RHD recording
% system. 
%
% inputs: 
%   basepath: path to SNLab data directory for a given session. 
% 
% optional:
%   video_idx: digitalin channel with video timestamps (default: 1 or intan
%   0)
%   event_idx: ditigalin channel where event timestamps will be saved (default is channel 2).
%   channel: channel where event timestamps are currently located (default
%   is channel 2 or intan channel 1).
%   amp_n_channels: number of channels in amplifier.dat (default 64)
%   aux_n_channels: number of channels in auxiliary.dat (default 3)
%   noise_epoch: matrix of start/end epochs for denoising (NOT CURRENTLY
%   WORKING LB 3/22). 
% output: 
%   - creates amplifier.dat and auxiliary.dat based on first and last timestamp in digtigalin.events.mat 
%     saves back to basepath.
%   - updates digitalin.events.mat to reflect trimmed amplifier/auxiliary
%     files
%   - renames digitalin, amplifier, and auxiliary dat files with '_old' to

% To-DO: 
% Add utility to replace epochs of dat with zeros (unplugs, etc). 
% Add option to trim based on epoch input. 
% Trim digitalin.dat, and time.dat. 

% LBerkowitz 03/22 
%
p = inputParser;
addParameter(p,'video_idx',1,@isnumeric) % default for processed data
addParameter(p,'event_idx',2,@isnumeric) % default for processed data
addParameter(p,'channel',2,@isnumeric) % default for processed data
addParameter(p,'amp_n_channels',64,@isnumeric) 
addParameter(p,'aux_n_channels',3,@isnumeric) 
addParameter(p,'noise_epochs',[],@isnumeric) 

parse(p,varargin{:});
video_idx = p.Results.video_idx;
event_idx = p.Results.event_idx;
channel = p.Results.channel;
amp_n_channels = p.Results.amp_n_channels;
aux_n_channels = p.Results.aux_n_channels;
noise_epochs = p.Results.noise_epochs; 


% load rhd to obtain recording parameters 
[~, ~, ~, ~,...
    ~, ~, frequency_parameters,~ ] = ...
    read_Intan_RHD2000_file_snlab(basepath);

% load aux and amp
amp = load_amplifier(basepath,amp_n_channels);
aux = load_auxiliary(basepath,aux_n_channels);
digitalIn = process_digitalin(basepath,frequency_parameters.board_dig_in_sample_rate);

% Pull start and end from digitalin channels that correspond to dat
%timestamp multiplied by sample rate
trim_dat_epoch(1) = digitalIn.timestampsOn{1, channel}(1)*...
    (frequency_parameters.board_dig_in_sample_rate);

trim_dat_epoch(2) = digitalIn.timestampsOff{1,channel}(end)*...
    (frequency_parameters.board_dig_in_sample_rate);

session_duration = (digitalIn.timestampsOff{1,channel}(end)...
    - digitalIn.timestampsOn{1, channel}(1))/60;

disp(['recording was ',num2str(session_duration),' minutes'])

% create batches 
calc_samples = trim_dat_epoch(2) - trim_dat_epoch(1);
batch = round(linspace(trim_dat_epoch(1),trim_dat_epoch(2),1000));

% write the data
write_amp(amp,basepath,batch)
write_aux(aux,basepath,batch)

% check size of amp compared to calcuated bytes
t = dir([basepath,filesep,'amplifier_.dat']);
samples = t.bytes/(amp_n_channels * 2); %int16 = 2 bytes
if samples == calc_samples
    disp([basepath,filesep,'amplifier_.dat ','size verified'])
else
    disp([basepath,filesep,'amplifier_.dat ','has different size than epoch'])
end


% check size of aux compared to calcuated bytes
t = dir([basepath,filesep,'auxiliary_.dat']);
samples = t.bytes/(aux_n_channels * 2); %int16 = 2 bytes
if samples == calc_samples
    disp([basepath,filesep,'auxiliary_.dat ','size verified'])
else
    disp([basepath,filesep,'auxiliary_.dat ','has different size than epoch'])
end

% adjust digitalin so time starts at zero
start_idx = (trim_dat_epoch(1)/frequency_parameters.board_dig_in_sample_rate);

digitalIn.timestampsOn{video_idx} = digitalIn.timestampsOn{1, video_idx} - start_idx;
digitalIn.timestampsOff{video_idx} = digitalIn.timestampsOff{1, video_idx} - start_idx;

digitalIn.timestampsOn{event_idx} = digitalIn.timestampsOn{1, channel} - start_idx;
digitalIn.timestampsOff{event_idx} = digitalIn.timestampsOff{1, channel} - start_idx;

% and recalculate other variables
for ii = 1:size(digitalIn.timestampsOn,2)
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

save([basepath,filesep,'digitalIn.events.mat'],'digitalIn');

% rename original files
system("rename " + basepath + filesep + "amplifier.dat" + " " + "amplifier_old.dat")
system("rename " + basepath + filesep +  "digitalin.dat" + " " + "digitalin_old.dat")
system("rename " + basepath + filesep +  "auxiliary.dat" + " " + "auxiliary_old.dat")

% rename trimmed files 
system("rename " + basepath + filesep + "amplifier_.dat" + " " + "amplifier.dat")
system("rename " + basepath + filesep +  "auxiliary_.dat" + " " + "auxiliary.dat")

end

function write_amp(amp,basepath,batch)
% skip if file is already created
if isfile([basepath,filesep,'amplifier_.dat'])
    disp([basepath,filesep,'amplifier_.dat ','already split'])
    return
end
amp_file = fopen(fullfile(basepath,'amplifier_.dat'),'w');
% loop though batches
for ii = 1:length(batch)-1
    disp(['batch ',num2str(batch(ii)+1),' to ',num2str(batch(ii+1)),...
        '   ',num2str(ii),' of ',num2str(length(batch)-1)])
    % write to disk
    fwrite(amp_file,amp.Data.mapped(:,batch(ii)+1:batch(ii+1)), 'int16');
end
fclose(amp_file);

end

function write_aux(aux,basepath,batch)
% LB 03/22

% skip if file is already created
if isfile([basepath,filesep,'auxiliary_.dat'])
    disp([basepath,filesep,'auxiliary_.dat ','already split'])
    return
end

% initiate file
aux_file = fopen([basepath,filesep,'auxiliary_.dat'],'w');

% loop though batches
for ii = 1:length(batch)-1
    disp(['batch ',num2str(batch(ii)+1),' to ',num2str(batch(ii+1)),...
        '   ',num2str(ii),' of ',num2str(length(batch)-1)])
    % write to disk
    fwrite(aux_file,aux.Data.mapped(:,batch(ii)+1:batch(ii+1)), 'uint16');
end
fclose(aux_file);
end