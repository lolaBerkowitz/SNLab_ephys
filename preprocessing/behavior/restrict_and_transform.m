function restrict_and_transform(basepath,varargin)
% restrict_and_transform takes xy coordinates from animal behavior file and
% transforms them to cm. 
%
% Assumes general behavior file and *maze_coords.csv is in basepath. 
p = inputParser;
p.addParameter('maze_size',[],@isnumeric)
p.addParameter('overwrite',false,@islogical)

p.parse(varargin{:});
maze_size = p.Results.maze_size;
overwrite = p.Results.overwrite;

% session basename to load animal behavior file and sessions file 
basename = basenameFromBasepath(basepath);

% load behavior file
load(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')
% check if units are in cm and if so return
if contains(behavior.position.units,'cm') && ~overwrite
    disp('Coordinates already scaled')
    return 
end

% Scales xy coordinates in behavior and saves back to behavior structure

% load sessions file
load(fullfile(basepath,[basename,'.session.mat']),'session')

% scales coordinates determined by known maze size and measurements in pixels
% from image
tracking.scale_coords(session,behavior,basepath);
% reload behavior file
load(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')

% restrict the coordinates and save back to behavior file
tracking.restrict(session,behavior,basepath);
end
 