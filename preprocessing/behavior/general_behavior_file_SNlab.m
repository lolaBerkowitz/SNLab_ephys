function general_behavior_file_SNlab(varargin)
% adapted from AYA lab general_behavior_file.m LB 03/22

% Outputs
p=inputParser;
addParameter(p,'basepath',pwd); % single or many basepaths in cell array or uses pwd
addParameter(p,'force_overwrite',false); % overwrite previously saved data (will remove custom fields)
addParameter(p,'save_mat',true); % save animal.behavior.mat
addParameter(p,'primary_coords_dlc',1:2); % deeplabcut tracking point to extract (extracts all, but main x and y will be this)
addParameter(p,'likelihood_dlc',.95); % deeplabcut likelihood threshold
addParameter(p,'smooth_factor',.1); % time in seconds to smooth over (default .1 or 100ms)
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

try 
    load(fullfile(basepath,'digitalIn.events.mat'))
catch
    disp('no digitalin.events found. Can''t run general behavior file for ephys')
    return
end

% check if file was already made
if force_overwrite
    disp('Overwriting previous runs')
elseif ~isempty(dir(fullfile(basepath,[basename,'.animal.behavior.mat'])))
    load([basepath,filesep,[basename,'.animal.behavior.mat']],'behavior');
    disp([basepath,filesep,[basename,'.animal.behavior.mat already detected. Loading file...']]);
    return
end

% run update_behavioralTracking to make sure dlc files and associated
% epochs are indicated in basename.session.behavioralTracking 
update_behavioralTracking('basepath',basepath)

% call extract_tracking which contains many extraction methods
[t,x,y,v,a,angle,units,source,fs,notes,extra_points,vidnames] = ...
    tracking.extract_tracking(basepath,primary_coords_dlc,likelihood_dlc,smooth_factor);

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
behavior.angle = angle';
behavior.trials = [];
behavior.trialsID = [];
behavior.states = [];
behavior.stateNames = [];
behavior.notes = notes;
behavior.epochs = session.epochs;
behavior.processinginfo.date = date;
behavior.processinginfo.vidnames = vidnames;
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


