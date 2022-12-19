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

% load behavior file
load(fullfile(basepath,[basename,'.animal.behavior.mat']))
% check if units are in cm and if so return
if ismember(behavior.position.units,'cm')
    disp('Coordinates already scaled')
    return 
end

% Scales xy coordinates in behavior and saves back to behavior structure

% load sessions file
load(fullfile(basepath,[basename,'.session.mat']))

% restrict the coordinates and save back to behavior file
tracking.restrict(session,behavior,basepath);

% scales coordinates determined by known maze size and measurements in pixels
% from image
tracking.scale_coords(session,behavior,basepath);

% save updates

end
 