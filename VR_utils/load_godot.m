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
parse(p,varargin{:});

savefile = p.Results.savefile;

% look for processed file and return if present
if ~isempty(dir(fullfile(basepath,'*vr_godot.csv')))
    disp('godot files processed. Loading csv')
    temp = dir(fullfile(basepath,'*vr_godot.csv'));
    vr_pos = readtable(fullfile(basepath,temp.name));
    return
end
% get file info
param = dir(fullfile(basepath,'*godotlogs.txt'));
vr_pos = table;

% run if file is found in basepath
if ~isempty(param)
    
    for lap = 1:length(param)
        temp_pos = readtable(fullfile(basepath,param(lap).name),'FileType','text');
        % identify trial from file name 
        temp_pos.lap_n = repmat(str2double(extractBetween(param(lap).name,'rep','_')),size(temp_pos,1),1);
        temp_pos.reward_n = repmat(str2double(extractBetween(param(lap).name,'reward','_')),size(temp_pos,1),1);
        temp_pos.mouse_id = repmat(str2double(extractBetween(param(lap).name,'mouse','_')),size(temp_pos,1),1);
        vr_pos = [vr_pos; temp_pos];
    end
    session_date = extractBefore(param(lap).name,'_linear');
    % add column labels
    vr_pos.Properties.VariableNames = {'yaw','pitch','roll','x',...
            'z','xz_angle','reward','lick','experiment_ts','lap_n','reward_n','mouse_id'};
    vr_pos = sortrows(vr_pos,'lap_n'); 
    if savefile
        writetable(vr_pos,fullfile(basepath,[session_date,'_vr_godot.csv']))
    end
else
    disp('No godotlogs found. Check basepath')
end 

end