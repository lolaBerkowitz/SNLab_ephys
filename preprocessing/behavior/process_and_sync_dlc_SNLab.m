function [tracking_struct,field_names] = process_and_sync_dlc_SNLab(varargin)
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
p.addParameter('primary_coords',1:2,@isnumeric); % which tracker point do you want
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
    load(fullfile(basepath,'digitalin.events.mat'));
    
    if exist('parsed_digitalIn','var')
        digitalIn = parsed_digitalIn;
    end
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
    
    [dlc_files,vidnames,ep_idx,frame_rate] =  tracking.get_tracking_info_from_session(session,basename);
    
else % create it and reload session and pull dlc,video, and epoch info
    % update session with dlc/video info
    update_behavioralTracking('basepath',basepath)
    session = loadSession(basepath,basename); % reload updated session file
    % pull tracking info
    [dlc_files,vidnames,ep_idx,frame_rate] =  tracking.get_tracking_info_from_session(session,basename);
    
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
    tempTracking{ii} = tracking.sync_ttl(basepath,video_ttl{ii},x,y,ts,fs,pulses_delta_range);
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
tracking_struct.position.x = x;
tracking_struct.position.y = y;
tracking_struct.folders = folder;
tracking_struct.samplingRate = samplingRate;
tracking_struct.timestamps = timestamps;
tracking_struct.description = description;
tracking_struct.vidnames = vidnames;

end




