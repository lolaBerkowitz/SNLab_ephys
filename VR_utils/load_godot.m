function vr_pos = load_godot(basepath,varargin)
% load_gadot
% 
% Input:
% basepath: A string representing the path to the directory containing 
%the GADOT experiment files. This directory should contain at least one .txt 
%file with the experiment data.
% 
% Output:
% A Matlab table containing the following columns: 
%'yaw', 'pitch', 'roll', 'x', 'z', 'xz_angle', 'reward', 'lick', 'experiment_ts'
% 
% Function Description:
% The load_gadot function reads in the data from the GADOT experiment 
% files located in the directory specified by basepath, 
% and returns a Matlab table containing the data.
% 
% The function assumes that the directory specified by basepath contains at least one .csv file with the experiment data. Each row of the data file should contain the following fields, in order:
% 
% yaw: a number representing the yaw angle of the head
% pitch: a number representing the pitch angle of the head
% roll: a number representing the roll angle of the head
% x: a number representing the x position of the head
% z: a number representing the z position of the head 
% xz_angle: a number representing the angle between the x-z plane and the head
% reward: a number representing whether a reward was given (1) or not (0)
% lick: a number representing whether a lick was detected (1) or not (0)
% experiment_ts: a number representing the timestamp of the experiment in milliseconds
% The function reads in the data from all the .csv files in the directory specified by basepath and combines them into a single Matlab table. If there are multiple files, the rows will be ordered by the experiment_ts column.
% 
% Example Usage:
% table_data = load_gadot('/path/to/gadot/files/')
% disp(table_data)

% Laura Berkowitz 2023

p = inputParser;
addParameter(p,'savefile',true,@islogical);
addParameter(p,'tags',{'lineartrack','openfield'},@iscell);

parse(p,varargin{:});

savefile = p.Results.savefile;
tags = p.Results.tags; 

% look for processed file and return if present
if ~isempty(dir(fullfile(basepath,'*vr_godot.csv')))
    disp('godot files processed. Loading csv')
    temp = dir(fullfile(basepath,'*vr_godot.csv'));
    vr_pos = readtable(fullfile(basepath,temp.name));
    return
end
% get file info
param = dir(fullfile(basepath,'*godotlogs.txt'));
if isempty(param)
    disp('No godot logs found.')
    vr_pos = table; 
    return
end

%% run if file is found in basepath
[param_new,task_type] = parse_param(param, tags);

for task = 1:length(param_new)
    param = param_new{task};
    task_name = task_type{task};
    
    [vr_pos,session_date] = load_text_files(param,basepath);

if savefile
    writetable(vr_pos,fullfile(basepath,[session_date,'_',task_name,'_vr_godot.csv']))
end

end


end

function [vr_pos,session_date] = load_text_files(param,basepath)
vr_pos = table;
for lap = 1:length(param)
    filename = param(lap).name;
    temp_pos = readtable(fullfile(basepath,filename),'FileType','text');
    % identify trial from file name
    temp_pos.lap_n = repmat(str2double(extractBetween(filename,'rep','_')),size(temp_pos,1),1);
    godot_date = regexp(filename, '(\d+)_(\d+)_(\d+)_(\d+)_(\d+)_(\d+)', 'match');
    temp_pos.session_date = repmat(godot_date,size(temp_pos,1),1);
%     temp_pos.reward_n = repmat(str2double(extractBetween(param(lap).name,'reward','_')),size(temp_pos,1),1);
%     temp_pos.mouse_id = repmat(str2double(extractBetween(param(lap).name,'mouse','_')),size(temp_pos,1),1);
    vr_pos = [vr_pos; temp_pos];
end

% add headers for each task type 
if contains(filename,'open')
    vr_pos.Properties.VariableNames = {'yaw','pitch','roll','x',...
        'z','xz_angle','reward','reward_x','reward_z','lick','experiment_ts','lap_n','godot_date'};

elseif contains(filename,'linear')
    % add column labels
    vr_pos.Properties.VariableNames = {'yaw','pitch','roll','x',...
        'z','xz_angle','reward','lick','experiment_ts','lap_n','godot_date'};
end

% sort rows in data. 
vr_pos = sortrows(vr_pos,{'godot_date','lap_n'});

% convert experiment timestamps to seconds
vr_pos.experiment_ts = vr_pos.experiment_ts/1000; % divide by 1000 as godot ts are in ms
session_date = vr_pos.godot_date{1};
end


function [param_new,task_type] = parse_param(param, tags)

% For different VR tasks, parse files into individual tasks 
for n = 1:length(tags)
    if any(contains({param.name},tags{n}))
        param_new{n} = param(contains({param.name},tags{n}));
        task_type{n} = tags{n};
    end
end

% remove empty cells 
param_new = param_new(~cellfun('isempty',param_new));
task_type = task_type(~cellfun('isempty',task_type));
end