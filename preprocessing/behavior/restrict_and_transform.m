function restrict_and_transform(basepath)
% restrict_and_transform takes xy coordinates from animal behavior file and
% transforms them to cm. 
%
% Assumes general behavior file and *maze_coords.csv is in basepath. 

% session basename to load animal behavior file and sessions file 
basename = basenameFromBasepath(basepath);

% load behavior file
load(fullfile(basepath,[basename,'.animal.behavior.mat']))

% load sessions file
load(fullfile(basepath,[basename,'.animal.behavior.mat']))


%% process coordinates in animal behavior file
start = [];
stop = [];
maze_size = [];

% loop through epochs to retrieve start/end used in restrictxy below
for ep = 1:length(session.behavioralTracking)
    epoch = session.behavioralTracking{1,ep}.epoch;
    start = [start,session.epochs{epoch}.startTime];
    stop = [stop,session.epochs{epoch}.stopTime];
    
    % load maze coords
    maze_coords_df = readtable(fullfile(basepath,...
        [extractBefore(session.behavioralTracking{1,ep}.notes,'.avi'),'_maze_coords.csv']));
    x_max = [x_max; max(maze_coords_df.x(ismember(maze_coords_df.object,'corner')))];
    x_min = [x_min; min(maze_coords_df.x(ismember(maze_coords_df.object,'corner')))];
    
    % gather maze size from behavior.epochs.maze_size 
    maze_size = [maze_size; behavior.epochs{1, epoch}.maze_size];

end

% gets index to restrict xy
good_idx = manual_trackerjumps(behavior.timestamps,...
    behavior.position.x,...
    behavior.position.y,...
    start,...
    stop,...
    basepath,'darkmode',false);

% for primary xy coordinates, restrict to boundaries provided by good_idx
behavior.position.x(~good_idx) = NaN;
behavior.position.y(~good_idx) = NaN;

% rescale coordinates
if lenth(x_max) > 1
    % two mazes average over scale factor
    for ep = 1:length(x_max)
        scale_factor(ep) = (x_max(ep) - x_min(ep))/maze_size(ep); %pixels/cm
    end
    % TEST HOW THIS WORKS WITH DIFFERENT MAZE SIZES SHOUL BE EQUIVALENT
    scale_factor = mean(scale_factor);
    
else
    scale_factor = (x_max(1) - x_min(1))/maze_size;
end

behavior.position.scale_factor = scale_factor;
coord_names = fieldnames(behavior.position);
% BREAK BELOW TO CHECK IF START/STOP CAN BE USED TO SCALE COODINATES 
% we'll only scale x and y for now, but will sor
for i = find(contains(coord_names,{'x','y'}))'
    behavior.position.(coord_names{i}) = behavior.position.(coord_names{i})/scale_factor;
end

behavior.position.units = 'cm';

save(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')
end