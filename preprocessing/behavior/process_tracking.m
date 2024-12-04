function process_tracking(basepath,varargin)
% processes dlc output and loads into animal behavior file for SNLab ephys
%
% Pipeline --> 
% 1. general_behavior_file creation -> 
% 2. trials_from_metadata -> 
% 3. maze_coords -> 
% 4. Restrictxy and scale coords -> 
% 5. saves have to general behavior file

% Dependencies: general_behavior_file_SNlab, update_behavior_from_metadata,
% get_maze_XY, restrict_and_transform
% 
% Assumptions: basepath contains videos of behavior, dlc output for each
% video, timestamps collected in digitalin.events{1,1}, and behavior
% session start/stop index in session.epochs. 
%
% Behavior metadata should be stored in metadata csv. 

% To-do: generate xy for trials as option, create behavior summary figure

% LB 2022
p = inputParser;
p.addParameter('overwrite_behavior',true,@islogical)
p.addParameter('config_path','C:\Users\schafferlab\github\SNLab_ephys\behavior\behavior_configs\')
p.addParameter('metadata_path','Y:\laura_berkowitz\behavior_metadata.csv')
p.addParameter('primary_coords_dlc',1:2,@isnumeric)


p.parse(varargin{:})
overwrite_behavior = p.Results.overwrite_behavior;
config_path = p.Results.config_path;
metadata_path = p.Results.metadata_path;
primary_coords = p.Results.primary_coords_dlc;

% skip if no tracking found
if isempty(dir([basepath,filesep,'*DLC*.csv'])) && isempty(dir([basepath,filesep,'*godot*.csv']))
    disp('No tracking found. Skipping session')
end

% skip if no digitalin events
try 
    load(fullfile(basepath,'digitalIn.events.mat'),'digitalIn')
catch
    disp('no digitalin.events found. Can''t sync tracking for general behavior file for ephys')
    return
end

% run main function
general_behavior_file_SNlab('basepath',basepath,'force_overwrite',overwrite_behavior,'primary_coords_dlc', primary_coords)

% update trials from metadata csv
update_behavior_from_metadata(metadata_path,'basepath',basepath)

% get maze coords
get_maze_XY('basepath',basepath,'config_path', config_path,'redo_rescale',true)

% sacle coordinates to cm 
tracking.scale_coords(basepath,overwrite_behavior);

% restrict coordintes to remove extramaze tracking points 
tracking.restrict(basepath,overwrite_behavior);
    
end



