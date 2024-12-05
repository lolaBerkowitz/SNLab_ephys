function process_behavior_folder(basepath,varargin)
% process_behavior_folder adds CellExplorer session, event, 
% and behavior files to folders that only contain behavior (no intan generated digitalin.events.

% Dependencies: 
    % SNLab_ephys/preprocessing/behavior
    % FFMPEG added as enviornmental variable so it can be called from CMD
    % Basepath added to a metadata.csv including columns basepath,
    %       trial_start_n, trial_stop_n, maze_width_cm, maze_length_cm,vidname
    % DLC coordinates (ear or electrode coordinates are 2,3)

% input: 
    % basepath: path to session folder including video and DLC files 
% output: 
    % .animal.behavior.mat - CellExplorer animal behavior file 
    % .restrictxy.mat - from neurocode/ephys_tools, provides boolean of
    %   excluded xy coordinates (outside of maze)
    % .session.mat - CellExplorer session file
    % *_maze_coords.csv - Coordinates of the maze boundaries in pixels and cm (after processing)
    %   for edge of maze.
    
% Example use case, 


    
%LB March 2024

% input parser
p = inputParser;
p.addParameter('metadata_path','Y:\laura_berkowitz\behavior_validation\appps1_cpp\metadata.csv',@ischar)
p.addParameter('overwrite',true,@islogical)
p.addParameter('redo_rescale',true,@islogical)


p.parse(varargin{:});
metadata_path = p.Results.metadata_path;
overwrite = p.Results.overwrite;
redo_rescale = p.Results.redo_rescale;

basename = basenameFromBasepath(basepath);

% skip if all files exist and overwrite is not equal to true
if exist([basepath,filesep,[basename,'.animal.behavior.mat']],'file')...
        & exist([basepath,filesep,[basename,'.restrictxy.mat']],'file')...
        & exist([basepath,filesep,[basename,'.session.mat']],'file')...
        & ~isempty(dir([basepath,filesep,'*_maze_coords.csv']))...
        & ~overwrite
    disp(['Folder: ', basepath, ' already processed.']) 
    return
end

% check is session file exists, if not make one 
if  ~exist([basepath,filesep,[basename,'.session.mat']],'file')
    session = sessionTemplate(basepath); 
    save(fullfile(basepath,[basename,'.session.mat']),'session')
end

% Make events from DLC
make_events('basepath',basepath);

% update epochs from digitalIn.events.mat
update_epochs('basepath',basepath,...
    'annotate',true,...
    'overwrite',overwrite,...
    'ttl_method',[])

% general behavior file
general_behavior_file_SNlab('basepath',basepath,'force_overwrite',overwrite,'smooth_factor',.2,'primary_coords_dlc',5);

% update behavior file from metadata csv
update_behavior_from_metadata(metadata_path,'basepath',basepath);

get_maze_XY('basepath',basepath,'config_path', 'C:\Users\schafferlab\github\SNLab_ephys\behavior\behavior_configs\','redo_rescale',redo_rescale,'overwrite',overwrite);

% sacle coordinates to cm 
tracking.scale_coords(basepath,overwrite);

% restrict coordintes to remove extramaze tracking points 
tracking.restrict(basepath,overwrite);
    

end