function [tracking,field_names] = process_and_sync_dlc_SNLab(varargin)
% Unpacks DLC and syncs with digitalin.events.mat timestamps
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
% L Berkowitz 03/2022 adapted some code by R Harvey

% to-do: 
% add functionality to choose primary coordinates 

% parse inputs
p = inputParser;
p.addParameter('basepath',pwd,@isfolder);
p.addParameter('video_channel_ttl',1,@isnumeric); % digitalin channel that contains video ttl
p.addParameter('primary_coords',1,@isnumeric); % which tracker point do you want
p.addParameter('likelihood',.95,@isnumeric); % tracking quality thres [0-1]
p.addParameter('pulses_delta_range',0.01,@isnumeric); % range for ttls

p.parse(varargin{:});
primary_coords = p.Results.primary_coords; % not used currently
basepath = p.Results.basepath;
video_channel_ttl = p.Results.video_channel_ttl;
likelihood = p.Results.likelihood;
pulses_delta_range = p.Results.pulses_delta_range;
basename = basenameFromBasepath(basepath);

% check for events, cannot process without them 
if exist(fullfile(basepath,'digitalin.events.mat'),'file') % only load if timestamps have been processed
    load(fullfile(basepath,'digitalin.events.mat'),'digitalIn');
else % if none, make them and update basename.session
    disp('processing events, one moment')
    % make digitalin.events.mat
    process_digitalin(basepath,session.extracellular.sr)
    % update epochs
    update_epochs('basepath',basepath,'annotate',true)
    % update behavioralTracking
    update_behavioralTracking('basepath',basepath)
end

% check for behavioralTracking field from session file and if present, grab
% tracking info (tracking file name, epoch index, and frame
% rate

% first load session file
session = loadSession(basepath,basename);

if isfield(session,'behavioralTracking')
    
    [dlc_files,vidnames,ep_idx,frame_rate] =  get_tracking_info_from_session(session,basename);
    
else % create it and reload session and pull dlc,video, and epoch info
    % update session with dlc/video info
    update_behavioralTracking('basepath',basepath)
    session = loadSession(basepath,basename); % reload updated session file
    % pull tracking info
    [dlc_files,vidnames,ep_idx,frame_rate] =  get_tracking_info_from_session(session,basename);
    
end

% grab video ttls
behav = 1;
for epoch = ep_idx
    % grab the video ttl timestamps within the epoch boundaries
    idx = digitalIn.timestampsOn{1, video_channel_ttl} >= session.epochs{1, epoch}.startTime  & ...
        digitalIn.timestampsOn{1, video_channel_ttl} <= session.epochs{1, epoch}.stopTime;
    video_ttl{behav} = digitalIn.timestampsOn{1, video_channel_ttl}(idx);
    behav = behav + 1;
    disp(['Epoch ', num2str(epoch), ' contains behavior tag. Grabbing video channel ttls for this epoch.'])
end

%% Sync video ttl with trackin file 
for ii = 1:length(dlc_files)
    disp(['processing file ',num2str(ii), ' of ',num2str(length(dlc_files))])
    
    % get frame rate of video
    fs = frame_rate(ii);
    
    % load csv with proper header
    opts = detectImportOptions(fullfile(basepath,dlc_files{ii}),'NumHeaderLines',2);
    df = readtable(fullfile(basepath,dlc_files{ii}),opts);
    
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
    x = df{:,x_col};
    y = df{:,y_col};
    
    % store tracking for each video file
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

% save data to ouput structure
tracking.position.x = x;
tracking.position.y = y;
tracking.folders = folder;
tracking.samplingRate = samplingRate;
tracking.timestamps = timestamps;
tracking.description = description;
tracking.vidnames = vidnames;

end

% Local functions 
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

function [tracking_files,vid_names,ep_idx,frame_rate] =  get_tracking_info_from_session(session,basename)
% extracts tracking information (tracking filename, video name, and epoch index
% basename.session.behavioralTracking

disp(['Checking for DLC files in ', basename, '.session.behavioralTracking'])

for i = 1:length(session.behavioralTracking)
    tracking_files{i} = session.behavioralTracking{1,i}.filenames;
    vid_names{i} = session.behavioralTracking{1,i}.notes;
    ep_idx(i) =  session.behavioralTracking{1,i}.epoch;
    frame_rate(i) = session.behavioralTracking{1,i}.framerate;
end

end


