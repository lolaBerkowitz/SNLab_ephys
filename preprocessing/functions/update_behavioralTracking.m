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
    'track','context','vr','pre_test','pairing_A','pairing_B'}); % save animal.behavior.mat
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

% check for videos in basepath
disp('Checking for behavior videos...')
vid_files = dir(fullfile(basepath,['*.','avi']));

% make sure earlier recorded videos are fist
[~,idx] = sort([vid_files.datenum]);
vid_files = vid_files(idx);

% check basepath for dlc tracking
dlc_files = get_dlc_files_in_basepath(basepath);
godot_files = dir(fullfile(basepath,'*vr_godot.csv'));

if isempty(vid_files) & isempty(godot_files)
    warning('No videos or godot logs found. behavioralTracking field not updated')
    return
elseif isempty(godot_files)
    warning('No VR godot logs found. behavioralTracking field not updated with godot')
elseif isempty(vid_files)
    warning('No video files found. behavioralTracking field not updated with videos.')
end

if isempty(dlc_files) & isempty(godot_files) & ~isempty(vid_files)
    warning('No DLC output found. Adding videos without tracking for now')
    dlc_files = vid_files;
end


beahve_files{1} = dlc_files;
beahve_files{2} = godot_files;


if ~isempty(dlc_files) | ~isempty(godot_files)
    
    session = add_behavior(session,beahve_files,vid_files,tags);
end

%     % update session with dlc and video files
%     session = add_videos(session,vid_files,dlc_files,tags);
% end
%
% if ~isempty(godot_files)
%     session = add_godot(session,godot_files,tags);
% end

% save data
save(fullfile(basepath,[basename, '.session.mat']),'session');

% verify videos added
if annotate
    disp('check epoch notes and verify videos were added to notes')
    gui_session(basepath)
end

end

function session = add_behavior(session,behave_files,vid_files,tags)

% if found update behavioralTracking field
disp('Tracking found!!')

% define basepath from basename.session
basepath = session.general.basePath;

% find epoch index with behavior tags
[ep_idx,epoch_name] = grab_epoch_index(session,tags);

% remove empty behave files 
 behave_files = behave_files(~cellfun('isempty',behave_files));
 behave_files = vertcat(behave_files{:});
% loop through videos
for i = 1:length(behave_files)
    
    % associate the epoch with the behavior file type     
    file = behave_files(i);
    filename = file(1).name;
        
    if contains(file(1).name,'_vr_godot')
        
        file_idx = find(contains({file.name},extractBefore(filename,'_vr_godot')));
        
        if length(file_idx) > 1
            warning('Multiple DLC outputs found for one video, keeping first index')
            file_idx = file_idx(1); % grab first one
        end
        
        session.behavioralTracking{1, i}.filenames = file(file_idx).name;
        session.behavioralTracking{1, i}.epoch  = ep_idx(find(contains(epoch_name,'vr')));
        session.behavioralTracking{1, i}.equipment  = 'VR_arduino_godot_controller_V3';
        session.behavioralTracking{1, i}.type  = 'MouseVR Godot Project V0.8';
        session.behavioralTracking{1, i}.notes  = filename;
        session.behavioralTracking{1, i}.framerate  = nan;
        session.behavioralTracking{1, i}.crop_params = [];
        
        % tracking is likely DLC
    else
        % this code associates a dlc file for a specific video, in case
        % there are multiple dlc/video files in basepath
        vidname = vid_files(i).name;
        dlc_idx = find(contains({vid_files.name},extractBefore({file.name},'DLC')));
        
        if length(dlc_idx) > 1
            warning('Multiple DLC outputs found for one video, keeping first index')
            dlc_idx = dlc_idx(1); % grab first one
        end
        %Pull up video
        videoObj = VideoReader(fullfile(basepath,vidname));
        
        % get model parameters for file
        crop_params = tracking.grab_dlc_crop(filename);
        
        % save to session file
        session.behavioralTracking{1, i}.filenames = filename;
        session.behavioralTracking{1, i}.epoch  = ep_idx(dlc_idx);
        session.behavioralTracking{1, i}.equipment  = 'FireFlyS';
        session.behavioralTracking{1, i}.type  = 'DLC';
        session.behavioralTracking{1, i}.notes  = vidname;
        session.behavioralTracking{1, i}.framerate  = videoObj.FrameRate;
        session.behavioralTracking{1, i}.crop_params = crop_params;
    end
    
    
end

end







function session = add_godot(session,godot_files,tags)

% if found update behavioralTracking field
disp('Godot tracking found!!')

% define basepath from basename.session
basepath = session.general.basePath;

% find epoch index with behavior tags
[ep_idx,epoch_name] = grab_epoch_index(session,tags);


ep_idx = ep_idx(contains(epoch_name,'vr'));
% loop through videos
for i = 1:length(godot_files)
    filename = godot_files(i).name;
    file_idx = find(contains({godot_files.name},extractBefore(filename,'_vr_godot')));
    
    if length(file_idx) > 1
        warning('Multiple DLC outputs found for one video, keeping first index')
        file_idx = file_idx(1); % grab first one
    end
    
    
    if contains(epoch_name(i),'vr')
        session.behavioralTracking{1, i}.filenames = godot_files(file_idx).name;
        session.behavioralTracking{1, i}.epoch  = ep_idx(i);
        session.behavioralTracking{1, i}.equipment  = 'VR_arduino_godot_controller_V3';
        session.behavioralTracking{1, i}.type  = 'MouseVR Godot Project V0.8';
        session.behavioralTracking{1, i}.notes  = filename;
        session.behavioralTracking{1, i}.framerate  = nan;
        session.behavioralTracking{1, i}.crop_params = [];
        
        % tracking is likely DLC
    else
        %Pull up video
        videoObj = VideoReader(fullfile(basepath,vidname));
        
        % get model parameters for file
        crop_params = tracking.grab_dlc_crop(dlc_files(dlc_idx).name);
        
        % save to session file
        session.behavioralTracking{1, i}.filenames = dlc_files(dlc_idx).name;
        session.behavioralTracking{1, i}.epoch  = ep_idx(i);
        session.behavioralTracking{1, i}.equipment  = 'FireFlyS';
        session.behavioralTracking{1, i}.type  = 'DLC';
        session.behavioralTracking{1, i}.notes  = vidname;
        session.behavioralTracking{1, i}.framerate  = videoObj.FrameRate;
        session.behavioralTracking{1, i}.crop_params = crop_params;
    end
    
    
end

end

function session = add_videos(session,vid_files,dlc_files,tags)

% if found update behavioralTracking field
disp('Videos found!!')

% define basepath from basename.session
basepath = session.general.basePath;

% find epoch index with behavior tags
[ep_idx,epoch_name] = grab_epoch_index(session,tags);

% for mupltiple videos
if length(vid_files) > 1
    %use timestamps appended to video name to infer epoch
    vid_ts = str2double(extractAfter(extractBefore({vid_files.name},'.avi'),'-'));
    
    % flip index if timestamp of video 2 is less then timestamp of
    % video 1
    if vid_ts(1) > vid_ts(2)
        ep_idx = flip(ep_idx);
    end
end

% loop through videos
for i = 1:length(vid_files)
    vidname = vid_files(i).name;
    dlc_idx = find(contains({dlc_files.name},extractBefore(vidname,'.avi')));
    
    if length(dlc_idx) > 1
        warning('Multiple DLC outputs found for one video, keeping first index')
        dlc_idx = dlc_idx(1); % grab first one
    end
    %Pull up video
    videoObj = VideoReader(fullfile(basepath,vidname));
    
    % get model parameters for file
    crop_params = tracking.grab_dlc_crop(dlc_files(dlc_idx).name);
    
    % save to session file
    session.behavioralTracking{1, i}.filenames = dlc_files(dlc_idx).name;
    session.behavioralTracking{1, i}.epoch  = ep_idx(i);
    session.behavioralTracking{1, i}.equipment  = 'FireFlyS';
    session.behavioralTracking{1, i}.type  = 'DLC';
    session.behavioralTracking{1, i}.notes  = vidname;
    session.behavioralTracking{1, i}.framerate  = videoObj.FrameRate;
    session.behavioralTracking{1, i}.crop_params = crop_params;
end


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
