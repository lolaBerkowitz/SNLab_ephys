function update_behavior_from_metadata(metadata_path,varargin)
% updates animal behavior file using info stored in metadata.csv

% input parser
p = inputParser;
addParameter(p,'basepath',pwd);

parse(p,varargin{:});
basepath = p.Results.basepath;

% read in metadata csv
df = readtable(metadata_path);
basename = basenameFromBasepath(basepath);
% skip if behavior file doesn't exist
try
    load(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')
catch
    return
end

% update trials in behavior file
update_trials(basepath,df)

% update maze size in session and in behavior file
update_maze_size(basepath,df)

% update pixel in session and in behavior file
update_pixel_distance(basepath,df)


end


function update_pixel_distance(basepath,df)
% updates behavior.trials from frames in df

% crate basename from basepath
basename = basenameFromBasepath(basepath);
% load session and animal behavior file
load(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior');

session = loadSession(basepath,basename);

% pull rows from df for this basepath
temp_df = df(contains(df.basename,basename),:);

% setup
vars = fieldnames(temp_df);
col_idx = contains(vars,{'pixel_distance','pixel_dist_reference'});

% loop through videos indicated in session.behavioralTracking
for ii = 1:length(session.behavioralTracking)
    
    % load epoch to get start/end for timestamps
    epoch = session.behavioralTracking{1,ii}.epoch;
    % grap video so we can find row index for a given video
    vidname = session.behavioralTracking{1,ii}.notes;
    row_idx = contains(temp_df.vidname,extractBefore(vidname,'.avi'));
        
    values = table2array(temp_df(row_idx,col_idx));
    behavior.epochs{1, epoch}.pixel_distance = values(1);
    behavior.epochs{1, epoch}.pixel_dist_reference = values(2);

    session.behavioralTracking{1,ii}.pixel_distance = values(1);
    session.behavioralTracking{1,ii}.pixel_dist_reference = values(2);


end
% save behavior file
save(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')
save(fullfile(basepath,[basename,'.session.mat']),'session')

end

function update_maze_size(basepath,df)
% updates behavior.trials from frames in df

% crate basename from basepath
basename = basenameFromBasepath(basepath);
% load session and animal behavior file
load(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior');

session = loadSession(basepath,basename);

% pull rows from df for this basepath
temp_df = df(contains(df.basename,basename),:);

% setup
vars = fieldnames(temp_df);
col_idx = contains(vars,{'maze_width_cm','maze_length_cm'});

% loop through videos indicated in session.behavioralTracking
for ii = 1:length(session.behavioralTracking)
    
    % load epoch to get start/end for timestamps
    epoch = session.behavioralTracking{1,ii}.epoch;
    % grap video so we can find row index for a given video
    vidname = session.behavioralTracking{1,ii}.notes;
    row_idx = contains(temp_df.vidname,extractBefore(vidname,'.avi'));
        
    session.behavioralTracking{1,ii}.maze_size = table2array(temp_df(row_idx,col_idx));

    behavior.epochs{1, epoch}.maze_size = table2array(temp_df(row_idx,col_idx));

end
% save behavior and session file
save(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')
save(fullfile(basepath,[basename,'.session.mat']),'session')
end

function update_trials(basepath,df)
% updates behavior.trials from frames in df

% crate basename from basepath
basename = basenameFromBasepath(basepath);
% load session and animal behavior file
load(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior');

session = loadSession(basepath,basename);

% pull rows from df for this basepath
temp_df = df(contains(df.basename,basename),:);

% setup
vars = fieldnames(temp_df);
col_idx = find(contains(vars,'trial'));
trial_ts = [];
trial_id = [];

t_n = 1;
% loop through videos indicated in session.behavioralTracking
for ii = 1:length(session.behavioralTracking)
    
    % exclude VR 
    if contains(session.behavioralTracking{1,ii}.filenames,'godot')
        continue
    end 
    % load epoch to get start/end for timestamps
    epoch = session.behavioralTracking{1,ii}.epoch;
    name = session.epochs{1,epoch}.name;
    % grap video so we can find row index for a given video
    vidname = session.behavioralTracking{1,ii}.notes;
    row_idx = contains(temp_df.vidname,extractBefore(vidname,'.avi'));
    
    % determine which columns have trial info
    frames = table2array(temp_df(row_idx,col_idx));
    
    % identify start and end frame columns
    start_idx = find(contains(vars(col_idx(~isnan(frames))),'start'));
    stop_idx = find(contains(vars(col_idx(~isnan(frames))),'stop'));
    
    % convert frames to seconds via behavior.timestamps for epoch
    start_ts = behavior.epochs{1, epoch}.startTime;
    stop_ts = behavior.epochs{1, epoch}.stopTime;
    
    % video ts from behavior file
    vid_ts = behavior.timestamps(behavior.timestamps >= start_ts & behavior.timestamps <= stop_ts);
    % frames are inicies for timestamps
    ts = vid_ts(frames(~isnan(frames)));
    for t_idx = 1:length(start_idx)
        temp_trial_ts(t_idx,:) = [ts(start_idx(t_idx)), ts(stop_idx(t_idx))];
    end
    trial_ts = [trial_ts; temp_trial_ts];
    clear temp_trial_ts
%     trial_name{t_n:t_n+length(ts(start_idx)')-1} = repmat(name,length(ts(start_idx)',1));
    for t = 1:length(ts(start_idx))
        trial_id{t_n} = [name,['_',num2str(t)]];
        t_n = t_n+1;
    end
    % update t_n so trials and trial_id will be same length

end

% Add to behavior file
behavior.trials = trial_ts;
behavior.trialsID = trial_id';
% save behavior file
save(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')
end
