function [tracking,field_names] = process_and_sync_dlc_SNLab(varargin)
% Unpacks DLC
%
% Run this after you have exported deeplabcut csv results, have
% basename.session.epochs labeled, and have digitalin.events.mat.
% 
% optional inputs: 
% basepath: location of SNLab session data (contains basename.session file,
%  digitalin.events.mat, and DLC output saved as CSV). 
% video_channel_ttl: digitalin channel that contains video ttl pulses (default channel 1). 
% primary_coords: cordinates of interest (default first dlc bodypart)
% 
%
%
% assumptions: 
% * Assumes video name has Spinview generated timestamps, and the dlc output
%   contains the video name (DLC default). This impacts the order of the video 
%   in the folder i.e. earlier videos first. Thus, the order of the DLC
%   output in the directory should correspond to the order of the ttl pulse
%   epoch (denoted by enviornment_name). 
%
%


% L Berkowitz 03/2022 adapted some code by R Harvey

p = inputParser;
p.addParameter('basepath',pwd,@isfolder);
p.addParameter('video_channel_ttl',1,@isnumeric);
p.addParameter('primary_coords',1,@isnumeric); % which tracker point do you want
p.addParameter('likelihood',.95,@isnumeric); % tracking quality thres [0-1]
p.addParameter('pulses_delta_range',0.01,@isnumeric); % range for ttls
p.addParameter('event_channel',2,@isnumeric); % default events epochs for SNlab is 2 for single session
p.addParameter('environment_name',{'linear_track','open_field','w_maze','y_maze','figure_8_maze'},@ischar);

p.parse(varargin{:});
primary_coords = p.Results.primary_coords; % not used currently
basepath = p.Results.basepath;
video_channel_ttl = p.Results.video_channel_ttl;
likelihood = p.Results.likelihood;
pulses_delta_range = p.Results.pulses_delta_range;
event_channel = p.Results.event_channel;
environment_name = p.Results.environment_name;

% load session
session = loadSession(basepath);

% check for events
if exist(fullfile(basepath,'digitalin.events.mat'),'file') % only load if timestamps have been processed
    load(fullfile(basepath,'digitalin.events.mat'));
    
    % for the sessions that had this saved incorrectly
    if exist('parsed_digitalIn','var')
        digitalIn = parsed_digitalIn;
        clear parsed_digitalIn
    end
    
else
    disp('processing events, one moment')
    % make digitalin.events.mat
    process_digitalin(basepath,session.extracellular.sr)
    
    % update epochs
    ii = 1;
    for i = 1:2:size(digitalin.timestampsOn{1, event_channel},1)-1 % by default 2nd column is events
        session.epochs{ii}.name =  char(i); % set label as empty
        session.epochs{ii}.startTime =  digitalin.timestampsOn{1, event_channel}(i);
        session.epochs{ii}.stopTime =  digitalin.timestampsOff{1, event_channel}(i+1);
        ii = ii+1;
    end
    disp('update session epochs with proper labels')
    gui_session(session)
end

% check for DLC csv
dlc_files = dir([basepath,filesep,'*DLC*.csv']);

if ~isempty(dlc_files)
    n_files = length(dlc_files);
    disp([num2str(n_files),' dlc files found.'])
    
    % grab video names from dlc file name
    vidnames = strcat(extractBefore({dlc_files.name},'DLC'),'.avi');
    
    disp('checking folder for videos')
    
    % find number of videos to check for dlc output
    vid_files = dir([basepath,filesep,'*.avi']);
    if sum(contains(vidnames,{vid_files.name})) == n_files
        disp('video files for dlc output found.')
    elseif sum(contains(vidnames,{vid_files.name}) )> 0
        disp('some files found, but not as many as dlc output. Check folder and add video for each dlc output')
        tracking = nan;
        field_names = nan;
        return
    else
        disp('No videos found. Check folder and add video for each dlc output')
        tracking = nan;
        field_names = nan;
        return
    end
end

%% find session epochs that contain behavior
behav = 1;
for epoch = 1:length(session.epochs)
    % grab the timestamps for the behavior epoch from the video channel
    if any(contains(fieldnames(session.epochs{1, epoch}),'environment'))
       if ismember(session.epochs{1, epoch}.environment,environment_name)
            idx = digitalIn.timestampsOn{1, video_channel_ttl} >= session.epochs{1, epoch}.startTime  & ...
                digitalIn.timestampsOn{1, video_channel_ttl} <= session.epochs{1, epoch}.stopTime;
            video_ttl{behav} = digitalIn.timestampsOn{1, video_channel_ttl}(idx);
            behav = behav + 1;
            disp(['Epoch ', num2str(epoch), ' contains behavior tag. Grabbing video channel ttls for this epoch.'])
       end
    else
        disp(['No behavior flag found in basename.session.epoch.evironment.',...
            'Enironment must contain flag indicated in the input enviornment_name'])
        return  
    end
end

%% Process dlc files below
for ii = 1:length({dlc_files.name})
    
    % load file
    file = dlc_files(ii).name;
    
    % load corresponding video
    video_file = vidnames(ismember(vidnames,strcat(extractBefore(file,'DLC'),'.avi')));
    obj = VideoReader(fullfile(basepath,video_file{1}));
    fs = obj.FrameRate;
    
    % load csv with proper header
    opts = detectImportOptions(fullfile(basepath,file),'NumHeaderLines',2);
    df = readtable(fullfile(basepath,file),opts);
    % get names of fields, these will be as long as tracker points
    % used times 3 because [x,y,likelihood]
    field_names = fields(df);
    % locate columns with [x,y,likelihood]
    x_col = find(contains(field_names,'x'));
    y_col = find(contains(field_names,'y'));
    likelihood_col = find(contains(field_names,'likelihood'));
    % filter out bad tracker points by likelihood thres
    for i = 1:length(x_col)
        idx = df{:,likelihood_col(i)} < likelihood;
        df{idx,x_col(i)} = NaN;
        df{idx,y_col(i)} = NaN;
    end
    ts = df{:,1}/fs;
    %             x = df{:,x_col(primary_coords)};
    %             y = df{:,y_col(primary_coords)};
    
    x = df{:,x_col};
    y = df{:,y_col};
    
    % find order of video in session (assumes trailing timestamps) 
    vid_files(ii)
    
    tempTracking{ii} = sync_ttl(basepath,video_ttl{ii},x,y,ts,fs,pulses_delta_range);
    trackFolder(ii) = ii;
end

% Concatenating tracking fields...
x = []; y = []; timestamps = []; folder = []; samplingRate = []; description = [];
for ii = 1:size(tempTracking,2)
    x = [x; tempTracking{ii}.position.x];
    y = [y; tempTracking{ii}.position.y];
    timestamps = [timestamps; tempTracking{ii}.timestamps];
    folder{ii} = tempTracking{ii}.folder;
    samplingRate = [samplingRate; tempTracking{ii}.samplingRate];
    description{ii} = tempTracking{ii}.description;
end

tracking.position.x = x;
tracking.position.y = y;
tracking.folders = folder;
tracking.samplingRate = samplingRate;
tracking.timestamps = timestamps;
end

function [tracking] = sync_ttl(basepath,video_ttl,x,y,ts,fs,pulses_delta_range)

% if ~exist(fullfile(folder,'digitalIn.events.mat'),'file')
%     digitalIn = getDigitalIn('all','folder',folder);
% end
% 
% load(fullfile(folder,'digitalIn.events.mat'))
% 
% Len = cellfun(@length, digitalIn.timestampsOn, 'UniformOutput', false);
% [~,idx] = max(cell2mat(Len));
% bazlerTtl = digitalIn.timestampsOn{idx};

%check for extra pulses of much shorter distance than they should
extra_pulses = diff(video_ttl)<((1/fs)-(1/fs)*pulses_delta_range);
video_ttl(extra_pulses) = [];

video_intan_diff = length(video_ttl) - size(x,1);

[x,y,ts,video_ttl] = match_video_frames_to_ttl(video_ttl,video_intan_diff,x,y,ts,fs);

[~,folder_name] = fileparts(basepath);
tracking.position.x = x;
tracking.position.y = y;
tracking.timestamps = video_ttl;
tracking.originalTimestamps = ts;
tracking.folder = folder_name;
tracking.samplingRate = fs;
tracking.description = '';
end

function [x,y,t,video_ttl] = match_video_frames_to_ttl(video_ttl,basler_intan_diff,x,y,t,fs)

% match basler frames con ttl pulses
if (length(video_ttl) == size(x,1)) || abs(basler_intan_diff)<=2 %assumes 1 frame could be cut at 0 and 1 frame at end
    disp('N of frames match!!');
elseif basler_intan_diff>0 && abs(basler_intan_diff)<fs
    disp([num2str(abs(length(video_ttl) - size(x,1))) ' of frames dont match, probably at the end of the recording']);
    video_ttl = video_ttl(1:size(x,1));
elseif basler_intan_diff<0 && abs(basler_intan_diff)<fs
    disp([num2str(abs(length(video_ttl) - size(x,1))) ' of frames dont match, probably at the beggining of the recording']);
    x = x(1:length(video_ttl),:);
    y = y(1:length(video_ttl),:);
elseif basler_intan_diff<0 && abs(basler_intan_diff)>fs
    disp([num2str(abs(size(x,1) - length(video_ttl)))...
        ' video frames without TTL... was the recording switched off before the camera? Cutting positions accordingly...']);
    x = x(1:length(video_ttl),:);
    y = y(1:length(video_ttl),:);
elseif abs(basler_intan_diff)>2*fs
    warning('More than 2 seconds missalignment in total in this session...will adjust to the closer one...');
    if basler_intan_diff>0
        video_ttl = video_ttl(1:size(x,1));
    else
        x = x(1:length(video_ttl),:);
        y = y(1:length(video_ttl),:);
    end
elseif isempty(video_ttl)
    video_ttl = t;
else
    warning('Unnoticed problem with Camera/Intan... I would go back and check both step by step');
end
end




