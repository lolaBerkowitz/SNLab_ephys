function general_behavior_file_SNlab(varargin)
% adapted from AYA lab general_behavior_file.m LB 03/22

% Outputs
p=inputParser;
addParameter(p,'basepath',pwd); % single or many basepaths in cell array or uses pwd
addParameter(p,'force_overwrite',false); % overwrite previously saved data (will remove custom fields)
addParameter(p,'save_mat',true); % save animal.behavior.mat
addParameter(p,'primary_coords_dlc',1:2); % deeplabcut tracking point to extract (extracts all, but main x and y will be this)
addParameter(p,'likelihood_dlc',.95); % deeplabcut likelihood threshold
addParameter(p,'smooth_factor',10); % n frames to smooth over (default 10 = 167ms for 60Hz)
% addParameter(p,'maze_size',30); % maze size in cm

parse(p,varargin{:});
basepath = p.Results.basepath;
force_overwrite = p.Results.force_overwrite;
save_mat = p.Results.save_mat;
primary_coords_dlc = p.Results.primary_coords_dlc;
likelihood_dlc = p.Results.likelihood_dlc;
smooth_factor = p.Results.smooth_factor;
% maze_size = p.Results.maze_size;

basename = basenameFromBasepath(basepath);

% check if file was already made
if exist([basepath,filesep,[basename,'.animal.behavior.mat']],'file') &&...
        ~force_overwrite
    load([basepath,filesep,[basename,'.animal.behavior.mat']],'behavior');
    disp([basepath,filesep,[basename,'.animal.behavior.mat already detected. Loading file...']]);
    return
end

% call extract_tracking which contains many extraction methods
[t,x,y,v,a,units,source,fs,notes,extra_points] = ...
    extract_tracking(basepath,primary_coords_dlc,likelihood_dlc,smooth_factor);

% load session file 
load([basepath,filesep,[basename,'.session.mat']]);

% package results
behavior.sr = fs;
behavior.timestamps = t';
behavior.position.x = x';
behavior.position.y = y';
behavior.position.z = []; % no z coords LB 06/22
behavior.position.linearized = [];
behavior.position.units = units;
behavior.speed = v';
behavior.acceleration = a';
behavior.trials = [];
behavior.states = [];
behavior.stateNames = [];
behavior.notes = notes;
behavior.epochs = session.epochs;
behavior.processinginfo.date = date;
behavior.processinginfo.function = 'general_behavioral_file_SNlab.mat';
behavior.processinginfo.source = source;

% deeplabcut will often have many tracking points, add them here
if ~isempty(extra_points)
    for field = fieldnames(extra_points)'
        field = field{1};
        behavior.position.(field) = extra_points.(field)';
    end
end

if save_mat
    save([basepath,filesep,[basename,'.animal.behavior.mat']],'behavior');
end
end

function [t,x,y,v,a,units,source,fs,notes,extra_points] = ...
    extract_tracking(basepath,primary_coords_dlc,likelihood_dlc,smooth_factor)

% initalize variables to pull
t = [];
x = [];
y = [];
v = [];
a = [];
units = [];
source = [];
notes = [];
extra_points = [];

% below are many methods on locating tracking data from many formats
[tracking,field_names] = process_and_sync_dlc_SNLab('basepath',basepath,...
    'primary_coords',primary_coords_dlc,...
    'likelihood',likelihood_dlc);

t = tracking.timestamps;
fs = 1/mode(diff(t));

x = tracking.position.x(:,primary_coords_dlc);
y = tracking.position.y(:,primary_coords_dlc);

if length(primary_coords_dlc) > 1
    % compute average point between two coords 
    x = median(x,2);
    y = median(y,2);
    [v, a,~] = linear_motion(x,y,fs,smooth_factor);
else
    [v, a,~] = linear_motion(x,y,fs,smooth_factor);
end

% multiple tracking points will likely exist, extract here
x_col = field_names(contains(field_names,'x'));
y_col = field_names(contains(field_names,'y'));
extra_points = struct();
for i = 1:length(x_col)
    extra_points.([x_col{i},'_point']) = tracking.position.x(:,i);
    extra_points.([y_col{i},'_point']) = tracking.position.y(:,i);
end

units = 'pixels';
source = 'deeplabcut';

if length(t) > length(x)
    t = t(1:length(x));
elseif length(x) > length(t)
    x = x(1:length(t));
    y = y(1:length(t));
    % adjust other tracker points
    for name = fields(extra_points)'
        extra_points.(name{1}) = extra_points.(name{1})(1:length(t));
    end
end

notes = ['primary_coords: ',num2str(primary_coords_dlc),...
    ', likelihood: ',num2str(likelihood_dlc)];

end
