function [tracking,field_names] = process_and_sync_dlc(varargin)
% Unpacks DLC
%
% Run this after you have exported deeplabcut csv results, have
% basename.session.epochs labeled, and have digitalin.events.mat.
%

% L Berkowitz 03/2022 adapted some code by RHarvey

p = inputParser;
p.addParameter('basepath',pwd,@isfolder);
p.addParameter('primary_coords',1,@isnumeric); % which tracker point do you want
p.addParameter('likelihood',.95,@isnumeric); % tracking quality thres [0-1]
p.addParameter('pulses_delta_range',0.01,@isnumeric); % range for ttls
p.addParameter('event_channel',2,@isnumeric); % default events epochs for SNlab is 2 for single session
p.addParameter('behavior_type',{'linear_track','open_field','w_maze','y_maze','figure_8_maze'},@ischar);

p.parse(varargin{:});
primary_coords = p.Results.primary_coords;
basepath = p.Results.basepath;
likelihood = p.Results.likelihood;
pulses_delta_range = p.Results.pulses_delta_range;
event_channel = p.Results.event_channel;
behavior_type = p.Results.behavior_type;

% load session
session = loadSession(basepath);

% grab basename
basename = basenameFromBasepath(basepath);

% check for events
if exist(fullfile(basepath,'digitalin.events.mat'),'file') % only load if timestamps have been processed
    load(fullfile(basepath,'digitalin.events.mat'));
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
        return
    else
        disp('No videos found. Check folder and add video for each dlc output')
    end
end

%% Process dlc files below

for ii = 1:length({dlc_files.name})
    file = dlc_files(ii).name;
    video_file = vidnames(ismember(vidnames,strcat(extractBefore(file,'DLC'),'.avi')));
%     obj = VideoReader(fullfile(basepath,video_file{1}));
%     fs = obj.FrameRate;
    
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
        
        tempTracking{ii} = sync_ttl(file.folder,x,y,ts,fs,pulses_delta_range);
        trackFolder(ii) = ii;
    end
end
end


% Concatenating tracking fields...
x = []; y = []; folder = []; samplingRate = []; description = [];
for ii = 1:size(tempTracking,2)
    x = [x; tempTracking{ii}.position.x];
    y = [y; tempTracking{ii}.position.y];
    folder{ii} = tempTracking{ii}.folder;
    samplingRate = [samplingRate; tempTracking{ii}.samplingRate];
    description{ii} = tempTracking{ii}.description;
end

tracking.position.x = x;
tracking.position.y = y;
tracking.folders = folder;
tracking.samplingRate = samplingRate;
tracking.timestamps = ts;
tracking.events.subSessions = subSessions;
tracking.events.subSessionsMask = maskSessions;
end

function [tracking] = sync_ttl(folder,x,y,ts,fs,pulses_delta_range)

if ~exist(fullfile(folder,'digitalIn.events.mat'),'file')
    digitalIn = getDigitalIn('all','folder',folder);
end
load(fullfile(folder,'digitalIn.events.mat'))

Len = cellfun(@length, digitalIn.timestampsOn, 'UniformOutput', false);
[~,idx] = max(cell2mat(Len));
bazlerTtl = digitalIn.timestampsOn{idx};

%check for extra pulses of much shorter distance than they should
extra_pulses = diff(bazlerTtl)<((1/fs)-(1/fs)*pulses_delta_range);
bazlerTtl(extra_pulses) = [];

basler_intan_diff = length(bazlerTtl) - size(x,1);

[x,y,ts,bazlerTtl] = match_basler_frames_to_ttl(bazlerTtl,basler_intan_diff,x,y,ts,fs);

[~,folder_name] = fileparts(folder);
tracking.position.x = x;
tracking.position.y = y;
tracking.timestamps = bazlerTtl;
tracking.originalTimestamps = ts;
tracking.folder = folder_name;
tracking.samplingRate = fs;
tracking.description = '';
end

function [x,y,t,bazlerTtl] = match_basler_frames_to_ttl(bazlerTtl,basler_intan_diff,x,y,t,fs)

% match basler frames con ttl pulses
if (length(bazlerTtl) == size(x,1)) || abs(basler_intan_diff)<=2 %assumes 1 frame could be cut at 0 and 1 frame at end
    disp('N of frames match!!');
elseif basler_intan_diff>0 && abs(basler_intan_diff)<fs
    disp([num2str(abs(length(bazlerTtl) - size(x,1))) ' of frames dont match, probably at the end of the recording']);
    bazlerTtl = bazlerTtl(1:size(x,1));
elseif basler_intan_diff<0 && abs(basler_intan_diff)<fs
    disp([num2str(abs(length(bazlerTtl) - size(x,1))) ' of frames dont match, probably at the beggining of the recording']);
    x = x(1:length(bazlerTtl),:);
    y = y(1:length(bazlerTtl),:);
elseif basler_intan_diff<0 && abs(basler_intan_diff)>fs
    disp([num2str(abs(size(x,1) - length(bazlerTtl)))...
        ' video frames without TTL... was the recording switched off before the camera? Cutting positions accordingly...']);
    x = x(1:length(bazlerTtl),:);
    y = y(1:length(bazlerTtl),:);
elseif abs(basler_intan_diff)>2*fs
    warning('More than 2 seconds missalignment in total in this session...will adjust to the closer one...');
    if basler_intan_diff>0
        bazlerTtl = bazlerTtl(1:size(x,1));
    else
        x = x(1:length(bazlerTtl),:);
        y = y(1:length(bazlerTtl),:);
    end
elseif isempty(bazlerTtl)
    bazlerTtl = t;
else
    warning('Unnoticed problem with Camera/Intan... I would go back and check both step by step');
end
end




