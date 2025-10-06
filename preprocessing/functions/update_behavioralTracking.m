function update_behavioralTracking(varargin)
% Adds behavioral files to session.beahviorTracking by looking for
% videos/dlc/godot output in session file
%
% Assumes basename.session.epochs.name have been annotated for each epoch
%
% inputs:
%   basepath: path where session data is kept. (Default is current working
%   directory)
%   annotate: boolean indicating whether user would like to check
%   basename.session via CellExplorer gui (Default is false)
%   tags: key words that finds behavioral epochs by matching
%   'basename.session.epochs.name'. (Recommend input as default is set for ease of
%   Berkowitz's app_ps1_ephys project {'base','learning','test','OF','open_field','hab','track','context'} )

% LB 2022
%input parser
p=inputParser;
addParameter(p,'basepath',pwd); % single or many basepaths in cell array or uses pwd
addParameter(p,'annotate',false); % save animal.behavior.mat
addParameter(p,'tags',{'base','learning','test','OF','open_field','morph',...
    'track','context','vr','pre_test','pairing_A','pairing_B','pairing',...
    'linear_track','post_test','y_maze','cheeseboard'}); % save animal.behavior.mat
addParameter(p,'force',false); % save animal.behavior.mat

parse(p,varargin{:});
basepath = p.Results.basepath;
annotate = p.Results.annotate;
tags = p.Results.tags;
force = p.Results.force;

basename = basenameFromBasepath(basepath);

% check if session contains behavioralTracking field
session = loadSession(basepath,basename);

if isfield(session,'behavioralTracking') && ~force
    disp('Session already contains field ''behavioralTracking''')
    return
end

% check for videos and godot files in basepath
disp('Checking for behavior videos...')
vid_files = dir(fullfile(basepath,['*.','avi']));
godot_files = dir(fullfile(basepath,'*vr_godot*.csv'));

% behavior files
behav_files = [vid_files; godot_files];

% order filenames by time of recording. Godot and videos have trailing
% timestamps.
behav_files = sort_by_trailing_ts({behav_files.name});

% Exit function if there are no behavioral files in the basepath
if isempty(behav_files)
    warning('No videos or godot logs found. behavioralTracking field not updated')
    return
end

% DLC files may accompany video files, so lets pull them from basepath
% check basepath for dlc tracking
dlc_files = get_dlc_files_in_basepath(basepath);
dlc_files = {dlc_files.name};

% if there are videos, but no dlc, we'll save vid_files in
% behavioralTracking field for now.
if isempty(dlc_files) && ~isempty(vid_files)
    warning('No DLC output found. Adding videos without tracking for now')
    dlc_files = {vid_files.name};
end

% update behavioralTracking field
session = update_field(session,behav_files,dlc_files,tags);

% save data
save(fullfile(basepath,[basename, '.session.mat']),'session');

% verify videos added
if annotate
    disp('check epoch notes and verify videos were added to notes')
    gui_session(basepath)
end

end

%% Local functions below

function session = update_field(session,behav_files,dlc_files,tags)


%% find epoch index with behavior tags
[ep_idx,~] = grab_epoch_index(session,tags);

% behavior files is sorted before input, so length of behavior files should
% equal the length of ep_idx. 
for i = 1:length(behav_files)
    if contains(behav_files(i),'godot')
        session = add_godot(session,behav_files(i),ep_idx(i),i);
    else
        session = add_dlc_files(session,dlc_files,behav_files(i),ep_idx(i),i); 
    end
    
end

end

function session = add_dlc_files(session,dlc_files,video_name,epoch,idx)

% if found update behavioralTracking field
disp('Adding DLC tracking now')

% Associate a video file for a specific video

if any(contains(dlc_files,'.avi'))
    dlc_file = video_name;
else
    dlc_file = dlc_files(contains(extractBefore(dlc_files,'DLC'),extractBefore(video_name,'.avi')));
    dlc_file = dlc_file{:};
end
%Pull up video
videoObj = VideoReader(fullfile(session.general.basePath,video_name{1}));

% get model parameters for file
crop_params = tracking.grab_dlc_crop(dlc_file);

% save to session file
session.behavioralTracking{1, idx}.filenames = dlc_file;
session.behavioralTracking{1, idx}.epoch  = epoch;
session.behavioralTracking{1, idx}.equipment  = 'FireFlyS';
session.behavioralTracking{1, idx}.type  = 'DLC';
session.behavioralTracking{1, idx}.notes  = video_name{1};
session.behavioralTracking{1, idx}.framerate  = videoObj.FrameRate;
session.behavioralTracking{1, idx}.crop_params = crop_params;


end

function session = add_godot(session,godot_files,epoch,idx)

% if found update behavioralTracking field
disp('Adding godot details now')

session.behavioralTracking{1, idx}.filenames = godot_files; % use first one
session.behavioralTracking{1, idx}.epoch  = epoch;
session.behavioralTracking{1, idx}.equipment  = 'VR_arduino_godot_controller_V3';
session.behavioralTracking{1, idx}.type  = 'MouseVR Godot Project V0.8';
session.behavioralTracking{1, idx}.notes  = godot_files;
session.behavioralTracking{1, idx}.framerate  = nan;
session.behavioralTracking{1, idx}.crop_params = [];

end

function [idx,epoch_name] = grab_epoch_index(session,tags)

for epoch = 1:length(session.epochs) % loop through epochs
    if any(contains(fieldnames(session.epochs{1, epoch}),'name'))
        epoch_name{epoch} = session.epochs{1, epoch}.name;
    else
        epoch_name{epoch} = 'none';
    end
    
end

idx = find(contains(epoch_name,tags,'IgnoreCase',true));
epoch_name = epoch_name(idx);
end

function filenames = sort_by_trailing_ts(filenames)
% takes in cell array of filenames that contain trailing timestamps ie
% filename-mmddyyyyhhmmss. Returns cell array sorted by earliest to latest


% Extract the timestamps from the filenames
if any(contains(filenames,'-0000'))
    timestamps = regexp(filenames, '-\d+-', 'match');
    timestamps = cellfun(@(x) extractAfter(x,'-'),timestamps,'UniformOutput',false);
    timestamps = cellfun(@(x) extractBefore(x,'-'),timestamps,'UniformOutput',false);
else 
    timestamps = regexp(filenames, '-\d+', 'match');
    timestamps = cellfun(@(x) extractAfter(x,'-'),timestamps,'UniformOutput',false);
end

timestamps = cellfun(@(x) str2double(x),timestamps,'UniformOutput',false);

% Sort the file list based on the timestamps
[~, sortedIdx] = sort( [timestamps{:}]);



filenames = filenames(sortedIdx);
end


