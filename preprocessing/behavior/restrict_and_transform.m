function restrict_and_transform(basepath,varargin)
% restrict_and_transform takes xy coordinates from animal behavior file and
% transforms them to cm. 
%
% Assumes general behavior file and *maze_coords.csv is in basepath. 
p = inputParser;
p.addParameter('maze_size',[],@isnumeric)

p.parse(varargin{:});
maze_size = p.Results.maze_size;

% session basename to load animal behavior file and sessions file 
basename = basenameFromBasepath(basepath);

% load sessions file
load(fullfile(basepath,[basename,'.session.mat']))

% restrict the coordinates and save back to behavior file
restrict(basepath,session);

% now scale the coordinates
[scale_factor_x,scale_factor_y] = scale_coords(session,basepath,maze_size);

% save maze coords as cm back to maze_coords.csv 
% loop through coords for each video
for ep = 1:length(session.behavioralTracking)
    maze_coords_df = readtable(fullfile(basepath,...
        [extractBefore(session.behavioralTracking{1,ep}.notes,'.avi'),'_maze_coords.csv']));
    maze_coords_df.x_cm = maze_coords_df.x/scale_factor_x;
    maze_coords_df.y_cm = maze_coords_df.y/scale_factor_y;
    
    writetable(maze_coords_df,fullfile(basepath,...
        [extractBefore(session.behavioralTracking{1,ep}.notes,'.avi'),'_maze_coords.csv']));
end

end

function restrict(basepath,session,behavior)

basename = basenameFromBasepath(basepath); 

% load behavior file
load(fullfile(basepath,[basename,'.animal.behavior.mat']))
start = [];
stop = [];
% gets index to restrict xy
for ep = 1:length(session.behavioralTracking)
    epoch = session.behavioralTracking{1,ep}.epoch;
    start = [start,session.epochs{epoch}.startTime];
    stop = [stop,session.epochs{epoch}.stopTime];
    
end

if ~isempty(dir(fullfile(basepath,[basename,'.restrictxy.mat'])))
    load(fullfile(basepath,[basename,'.restrictxy.mat']))
else
    good_idx = manual_trackerjumps(behavior.timestamps,...
        behavior.position.x,...
        behavior.position.y,...
        start,...
        stop,...
        basepath,'darkmode',false);
end
% for primary xy coordinates, restrict to boundaries provided by good_idx
behavior.position.x(~good_idx) = NaN;
behavior.position.y(~good_idx) = NaN;
save(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')

end

function [scale_factor_x,scale_factor_y] = scale_coords(session,basepath,maze_size)


basename = basenameFromBasepath(basepath); 

% load behavior file
load(fullfile(basepath,[basename,'.animal.behavior.mat']))
% Scales xy coordinates in behavior and saves back to behavior structure

if isempty(maze_size)
    maze_size = [];
end


x_max = [];
x_min = [];
y_max = [];
y_min = [];


% loop through epochs to retrieve start/end used in restrictxy below
for ep = 1:length(session.behavioralTracking)
    epoch = session.behavioralTracking{1,ep}.epoch;
    
    % load maze coords
    maze_coords_df = readtable(fullfile(basepath,...
        [extractBefore(session.behavioralTracking{1,ep}.notes,'.avi'),'_maze_coords.csv']));
    x_max = [x_max; max(maze_coords_df.x(ismember(maze_coords_df.object,'corner')))];
    x_min = [x_min; min(maze_coords_df.x(ismember(maze_coords_df.object,'corner')))];
    y_max = [y_max; max(maze_coords_df.y(ismember(maze_coords_df.object,'corner')))];
    y_min = [y_min; min(maze_coords_df.y(ismember(maze_coords_df.object,'corner')))];
    
    % gather maze size from behavior.epochs.maze_size 
    if ~isempty(behavior.epochs{1, epoch}.maze_size)
        maze_size = [maze_size; behavior.epochs{1, epoch}.maze_size];
    end

end

% rescale coordinates
if length(x_max) > 1
    % two mazes average over scale factor
    for ep = 1:length(x_max)
        scale_factor_x(ep) = (x_max(ep) - x_min(ep))/maze_size(ep); %pixels/cm
        scale_factor_y(ep) = (y_max(ep) - y_min(ep))/maze_size(ep); %pixels/cm
    end
    % TEST HOW THIS WORKS WITH DIFFERENT MAZE SIZES SHOUL BE EQUIVALENT
    scale_factor_x = mean(scale_factor_x);
    scale_factor_y = mean(scale_factor_y);
    
else
    scale_factor_x = (x_max(1) - x_min(1))/maze_size;
    scale_factor_y = (y_max(1) - y_min(1))/maze_size;

end

if ismember(behavior.position.units,'cm')
    disp('Coordinates already scaled')
    return 
end
    

coord_names = fieldnames(behavior.position);

% BREAK BELOW TO CHECK IF START/STOP CAN BE USED TO SCALE COODINATES 
% we'll only scale x and y for now, but will sor
for i = find(contains(coord_names,{'x'}))'
    behavior.position.(coord_names{i}) = behavior.position.(coord_names{i})/scale_factor_x;
end

for i = find(contains(coord_names,{'y'}))'
    behavior.position.(coord_names{i}) = behavior.position.(coord_names{i})/scale_factor_y;
end
behavior.position.units = 'cm';
save(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')

end