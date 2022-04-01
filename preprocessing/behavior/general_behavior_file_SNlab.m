function general_behavior_file_SNlab(varargin)
% adapted from AYA lab general_behavior_file.m
% Outputs
p=inputParser;
addParameter(p,'basepath',pwd); % single or many basepaths in cell array or uses pwd
addParameter(p,'force_overwrite',false); % overwrite previously saved data (will remove custom fields)
addParameter(p,'force_run',true); % run even if animal.behavior already exists
addParameter(p,'save_mat',true); % save animal.behavior.mat
addParameter(p,'primary_coords_dlc',1:2); % deeplabcut tracking point to extract (extracts all, but main x and y will be this)
addParameter(p,'likelihood_dlc',.95); % deeplabcut likelihood threshold

parse(p,varargin{:});
basepaths = p.Results.basepath;
fs = p.Results.fs;
force_overwrite = p.Results.force_overwrite;
force_run = p.Results.force_run;
save_mat = p.Results.save_mat;
primary_coords_dlc = p.Results.primary_coords_dlc;
likelihood_dlc = p.Results.likelihood_dlc;

% initalize variables to pull
t = [];
x = [];
y = [];
v = [];
trials = [];
units = [];
source = [];
notes = [];
extra_points = [];

extract_tracking(basepath,primary_coords_dlc,likelihood_dlc)


end

function [t,x,y,z,v,units,source,fs,notes,extra_points] = ...
    extract_tracking(basepath,primary_coords_dlc,likelihood_dlc)l

% initalize variables to pull
t = [];
x = [];
y = [];
v = [];
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

x = tracking.position.x(:,primary_coords);
y = tracking.position.y(:,primary_coords);

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

notes = ['primary_coords: ',num2str(primary_coords),...
    ', likelihood: ',num2str(likelihood_dlc)];
notes = {notes,tracking.notes};
end
