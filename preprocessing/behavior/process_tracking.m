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
p.addParameter('maze_coords',true,@islogical)
p.addParameter('config_path','C:\Users\schafferlab\github\SNLab_ephys\behavior\behavior_configs\')
p.addParameter('metadata_path','Y:\laura_berkowitz\app_ps1_ephys\behavior\object_location\object_location_metadata.csv')

p.parse(varargin{:})
maze_coords = p.Results.maze_coords;
config_path = p.Results.config_path;
metadata_path = p.Results.metadata_path;

basename = basenameFromBasepath(basepath);
session = loadSession(basepath,basename);

% run main function
general_behavior_file_SNlab('basepath',basepath)

% update trials from metadata csv
update_behavior_from_metadata(metadata_path,'basepath',basepath)

% get maze coords
if maze_coords
    config = get_behavior_config(session,config_path);
    get_maze_XY('basepath',basepath,'config_path', fullfile(config_path,config))
end

% restrict and transform primary coordinates to maze and convert to cm
restrict_and_transform(basepath)

end


