function general_behavior_file_SNlab(varargin)
% adapted from AYA lab general_behavior_file.m LB 03/22

% Outputs
p=inputParser;
addParameter(p,'basepath',pwd); % single or many basepaths in cell array or uses pwd
addParameter(p,'force_overwrite',false); % overwrite previously saved data (will remove custom fields)
addParameter(p,'save_mat',true); % save animal.behavior.mat
addParameter(p,'primary_coords_dlc',3); % deeplabcut tracking point to extract (extracts all, but main x and y will be this)
addParameter(p,'likelihood_dlc',.90); % deeplabcut likelihood threshold

parse(p,varargin{:});
basepath = p.Results.basepath;
force_overwrite = p.Results.force_overwrite;
save_mat = p.Results.save_mat;
primary_coords_dlc = p.Results.primary_coords_dlc;
likelihood_dlc = p.Results.likelihood_dlc;

basename = basenameFromBasepath(basepath);
session = loadSession(basepath,basename);

% check if file was already made
if force_overwrite
    disp('Overwriting previous runs')
elseif ~isempty(dir(fullfile(basepath,[basename,'.animal.behavior.mat'])))
    disp([basepath,filesep,[basename,'.animal.behavior.mat already detected. Loading file...']]);
    load([basepath,filesep,[basename,'.animal.behavior.mat']],'behavior');
    return
end

% run update_behavioralTracking to make sure dlc files and associated
% epochs are indicated in basename.session.behavioralTracking
update_behavioralTracking('basepath',basepath,'force',force_overwrite)
session = loadSession(basepath,basename);


% get tracking files
for i = 1:length(session.behavioralTracking)
    tracking_files{i} = session.behavioralTracking{1, i}.filenames;
end

% extract tracking for godot tracking
if any(contains({tracking_files{:}},'godot'))
    [t,x,y,v,a,angle,units,source,fs,notes,extra_points,vidnames] = ...
        tracking.extract_godot_tracking(basepath);
end


% extract tracking for Deeplabcut 
if any(contains({tracking_files{:}},'DLC'))
    [t_dlc,x_dlc,y_dlc,v_dlc,a_dlc,angle_dlc,units_dlc,source_dlc,fs_dlc,notes_dlc,extra_points_dlc,vidnames_dlc] = ...
        tracking.extract_tracking(basepath,primary_coords_dlc,likelihood_dlc);
    
    % deeplabcut will often have many tracking points, add them here
    if ~isempty(extra_points_dlc)
        for field = fieldnames(extra_points_dlc)'
            behavior.position.(field{1}) = extra_points_dlc.(field{1})';
        end
    end
    
end

% transpose timestamps
if exist('t','var') && isrow(t)
    t = t';
end

% transpose timestamps from DLC
if exist('t_dlc','var') && isrow(t_dlc)
    t_dlc = t_dlc';
end

% combine both
if exist('t','var') && exist('t_dlc','var')
    
    if min(t) < min(t_dlc)
        % concatenate
        t = [t; t_dlc];
        x = [x; x_dlc];
        y = [y; y_dlc];
        v = [v; v_dlc];
        a = [a; a_dlc];
        angle = [angle; angle_dlc];
        units = {units, units_dlc};
        source = {source; source_dlc};
        fs = {fs; fs_dlc};
        notes = {notes; notes_dlc};
        vidnames = {vidnames; vidnames_dlc};
        extra_points = {extra_points; extra_points_dlc};
        
    else
        % concatenate
        t = [t_dlc; t];
        x = [x_dlc; x];
        y = [y_dlc; y];
        v = [v_dlc; v];
        a = [a_dlc; a];
        angle = [angle_dlc; angle];
        units = {units_dlc;units};
        source = {source_dlc;source};
        fs = {fs_dlc;fs};
        notes = {notes_dlc;notes};
        vidnames = {vidnames_dlc;vidnames};
        extra_points = {extra_points_dlc;extra_points};
    end
end

% for DLC only 
if ~exist('t','var') && exist('t_dlc','var')
        % rename
    t = t_dlc;
    x = x_dlc;
    y = y_dlc;
    v = v_dlc;
    a = a_dlc;
    angle = angle_dlc;
    units = {units_dlc};
    source = {source_dlc};
    fs = fs_dlc;
    notes = {notes_dlc};
    vidnames = {vidnames_dlc};
    extra_points = {extra_points_dlc};
    
end

% load session file
load([basepath,filesep,[basename,'.session.mat']],'session');

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

if save_mat
    save([basepath,filesep,[basename,'.animal.behavior.mat']],'behavior');
end
end


