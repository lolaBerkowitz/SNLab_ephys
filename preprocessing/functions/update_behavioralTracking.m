function update_behavioralTracking(varargin)
% Adds video files to session.beahviorTracking by looking for videos/dlc
% output in session file
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
addParameter(p,'tags',{'base','learning','test','OF','open_field','morph','track','context'}); % save animal.behavior.mat

parse(p,varargin{:});
basepath = p.Results.basepath;
annotate = p.Results.annotate;
tags = p.Results.tags;
basename = basenameFromBasepath(basepath);

% check for videos in basepath
disp('Checking for behavior videos...')
vid_files = dir(fullfile(basepath,['*.','avi']));

% make sure earlier recorded videos are fist
[~,idx] = sort([vid_files.datenum]);
vid_files = vid_files(idx);

% check basepath for dlc tracking
dlc_files = dir([basepath,filesep,'*DLC*.csv']); % check for dlc output

session = loadSession(basepath,basename);

if isempty(vid_files)
    warning('No videos found. behavioralTracking field not updated. Exiting function')
    return
end

if isempty(dlc_files)
    warning('No DLC output found. Adding videos without tracking for now')
    dlc_files = vid_files;
end

% update session with dlc and video files 
session = main(session,vid_files,dlc_files,tags); 

% save data
save(fullfile(basepath,[basename, '.session.mat']),'session');

% verify videos added
if annotate
    disp('check epoch notes and verify videos were added to notes')
    gui_session('basepath',basepath)
end

end

function session = main(session,vid_files,dlc_files,tags)

% if found update behavioralTracking field
disp('Videos found!!')

% define basepath from basename.session
basepath = session.general.basePath;

% find epoch index with behavior tags
[ep_idx,epoch_name] = grab_epoch_index(session,tags);

% for mupltiple videos
if length(epoch_name) > 1
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
    
    % save to session file
    session.behavioralTracking{1, i}.filenames = dlc_files(dlc_idx).name;
    session.behavioralTracking{1, i}.epoch  = ep_idx(i);
    session.behavioralTracking{1, i}.equipment  = 'FireFlyS';
    session.behavioralTracking{1, i}.type  = 'DLC';
    session.behavioralTracking{1, i}.notes  = vidname;
    session.behavioralTracking{1, i}.framerate  = videoObj.FrameRate;
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
