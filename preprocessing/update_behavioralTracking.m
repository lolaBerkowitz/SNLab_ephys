function update_behavioralTracking(varargin)
% Adds video files to session.beahviorTracking

%input parser
p=inputParser;
addParameter(p,'basepath',pwd); % single or many basepaths in cell array or uses pwd
addParameter(p,'annotate',false); % save animal.behavior.mat
addParameter(p,'tags',{'base','learning','test','OF','open_field','hab','track','context'}); % save animal.behavior.mat

% addParameter(p,'maze_size',30); % maze size in cm

parse(p,varargin{:});
basepath = p.Results.basepath;
annotate = p.Results.annotate;
tags = p.Results.tags;

basename = basenameFromBasepath(basepath);
disp('Checking for behavior videos...')

% check basepath for dlc videos
session = loadSession(basepath,basename);

% check for videos in basepath
disp('Checking for behavior videos...')
vid_files = dir(fullfile(basepath,['*.','avi']));
dlc_files = dir([basepath,filesep,'*DLC*.csv']); % check for dlc output

% if found update behavioralTracking field
if ~isempty(vid_files)
    disp('Videos found')
    
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
       
        %Pull up video
        videoObj = VideoReader(fullfile(basepath,vidname));
        
        % save to session file
        session.behavioralTracking{1, dlc_idx}.filenames = dlc_files(dlc_idx).name;
        session.behavioralTracking{1, dlc_idx}.epoch  = ep_idx(i);
        session.behavioralTracking{1, dlc_idx}.equipment  = 'FireFlyS';
        session.behavioralTracking{1, dlc_idx}.type  = 'DLC';
        session.behavioralTracking{1, dlc_idx}.notes  = vidname;
        session.behavioralTracking{1, dlc_idx}.framerate  = videoObj.FrameRate;
    end
    
    save(fullfile(basepath,[basename, '.session.mat']),'session');
    % verify videos added
    if annotate
        disp('check epoch notes and verify videos were added to notes')
        gui_session(basepath)
    end
end

% save the file
save(fullfile(basepath,[basename, '.session.mat']),'session');


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
