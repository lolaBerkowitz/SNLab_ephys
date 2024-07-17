function update_behavior_from_metadata(metadata_path,varargin)
% updates animal behavior file using info stored in metadata.csv

% input parser
p = inputParser;
addParameter(p,'basepath',pwd);
addParameter(p,'batch',false); % flag for batch processing
addParameter(p,'data_folder',[]); % directory for batch

parse(p,varargin{:});
basepath = p.Results.basepath;
batch = p.Results.batch;
data_folder = p.Results.data_folder; % for batch in progress

% read in metadata csv
df = readtable(metadata_path);
basename = basenameFromBasepath(basepath);
% skip if behavior file doesn't exist
try
    load(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')
catch
    return
end

if batch
    disp('batch processing option a work in progress')
    return
    
end

% update trials
update_trials(basepath,df)

% update maze size
update_maze_size(basepath,df)


end



% Batch function in progress
function main_batch(data_folder,df)
% loop through and update animal behavior file for all subfolders in
% data_folder
%
basenames = unique(df.basename);
% loop through unique basenamesl
for i = 1:length(basenames)
    % obtain subject dir from basename
    basename = basenames{i};
    
    sub = split(basename,'_');
    sub = sub{1};
    
    % create basepath
    basepath = [data_folder,filesep,sub,filesep,basename];
    session = loadSession(basepath,basename);
    
    % load animal behavior file
    if exist(fullfile(basepath,[basename,'.animal.behavior.mat']),'file')
        load(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')
    else
        try
            general_behavior_file_SNlab('basepath',basepath)
            load(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')
        catch
            disp('Failed - check DLC outputs are present')
            continue
        end
        
    end
    
    if ~isempty(behavior.trials)
        continue
    else
    %% update behavior.trials from trial_start/stop columns
    update_trials(basepath,df)
    end
    % update behavior.epochs.maze_size from maze_size column
    update_maze_size(basepath,df)
    
    
end
end


function update_maze_size(basepath,df)
% updates behavior.trials from frames in df

% crate basename from basepath
basename = basenameFromBasepath(basepath);
% load session and animal behavior file
load(fullfile(basepath,[basename,'.animal.behavior.mat']));

session = loadSession(basepath,basename);

% pull rows from df for this basepath
temp_df = df(contains(df.basename,basename),:);

% setup
vars = fieldnames(temp_df);
col_idx = find(contains(vars,{'maze_width_cm','maze_length_cm'}));

% loop through videos indicated in session.behavioralTracking
for ii = 1:length(session.behavioralTracking)
    
    % load epoch to get start/end for timestamps
    epoch = session.behavioralTracking{1,ii}.epoch;
    % grap video so we can find row index for a given video
    vidname = session.behavioralTracking{1,ii}.notes;
    row_idx = contains(temp_df.vidname,extractBefore(vidname,'.avi'));
        
    behavior.epochs{1, epoch}.maze_size = table2array(temp_df(row_idx,col_idx));

end
% save behavior file
save(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')
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
