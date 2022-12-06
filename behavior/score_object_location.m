function score_object_location(varargin)
% Computes measures of object expoloration for object-location behavior
% sessions. 

% inputs
p = inputParser;
p.addParameter('basepath',pwd,@islogical)

p.parse(varargin{:})
basepath = p.Results.basepath;

% Load behavior file from basepath and check for (epochs/trials) 
% load maze_coords.csv and obtain xy position of objects 
load(fullfile(basepath,[basename,'.animal.behavior.mat']))
coord_paths = dir(fullfile(basepath,'*_maze_coords.csv'));

% check input dependencies 
if ~exist('behavior','var')
    error('No behavior file found. Check basepath for basename.animal.behavior.mat')
end 

if isempty(coord_paths)
    error('No maze_coords found. Run get_maze_XY.m')
end 


% create object boundary 

% determine moved object by measuring differences in epoch 1 vs epoch 2
% position. 

% find time within object boundary for each object as a function of time in
% session (5 minute bins for entire session)

% compute discrimination index for each bin 

% compute other measures of locomotion in each bin (speed, path length,
% time spent moving (movement greater than 3cm/s). 

% Save data to df 

end