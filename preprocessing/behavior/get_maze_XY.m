function get_maze_XY(varargin)
% getXY allows user to obtain xy coordinates of maze corners and objects
% (center/edge) from a video. Saves csv for each folder.
%
%   ASSUMPTIONS:
%    Videos must be in format compatible with Video Reader (see matlab
%   documentation for compatibilities given your computers OS). Recommend .AVI
%   as its compatible for most OS.

%
%   Subdirectories names become subID in saved table.
%   INPUTS (Optional) :
%       basepath: session folder containing video(s).
%       vid_type: file extension of video (default '.avi')
%       vid_time: time in seconds to load video.
%       ephys: boolean indicating ephys session. Assumes CellExplorer
%       processing has been done.
%       multi_chamber: boolean indicating if multiple behavior chambers are
%       in frame.
%   OUTPUT:
%       xycoords --> table containing max min values for each video file.
%
% examples:
%
% Behavior file with multi_chamber videos:
%   getXY('basepath',basepath,'multi_chamber',true)
%
% To-do:
%  add ability to scroll through frames for start/end trials
%
% Main dependencies: CellExplorer repo, general_behavior_file_SNLab, update_behavioralTracking
% LB 2022

% input parser
p = inputParser;
p.addParameter('basepath',pwd,@isfolder);
p.addParameter('vid_type','.avi',@ischar);
p.addParameter('vid_time',300,@isnumeric); % time of video to load in seconds
p.addParameter('config_path','C:\Users\schafferlab\github\SNLab_ephys\behavior\behavior_configs',@islogical); % time of video to load in seconds


p.parse(varargin{:})
basepath = p.Results.basepath; % not used currently
vid_type = p.Results.vid_type; % not used currently
vid_time = p.Results.vid_time; % not used currently
config_path = p.Results.config_path;

% find all video files indicated by vid_type
basename = basenameFromBasepath(basepath);
session = loadSession(basepath,basename);

if ~isfield(session,'behavioralTracking')
    warning('No tracking items found. No maze coords to return.')
    return
end

% check if behavior file exsists and if not make one
if ~exist([basepath,filesep,[basename,'.animal.behavior.mat']],'file')
    general_behavior_file_SNlab(basepath)
end

% get coords
main(session, config_path, vid_time)


end

function main(session, config_path, vid_time)
% runs main process of looping through videos in basepath, pulling up
% image via local grab_coords function which allows user to collect coordinate data
% and outputs coords for object A center, object A edge, object B center,
% object B edge, corner center A-D (inner corner of maze).
%

% LB 2022
basepath = session.general.basePath;
basename = basenameFromBasepath(basepath);

% load animal behavior file
load(fullfile(basepath,[basename,'.animal.behavior.mat']))


% loop through video
for file = 1:length(session.behavioralTracking) %loop through folders containing subject videos
    
    vid_path = fullfile(basepath,session.behavioralTracking{1,file}.notes);
    
    %Pull up video
    videoObj = VideoReader(vid_path,'CurrentTime',vid_time); % load video starting at vid_time
    
    epoch = session.behavioralTracking{1,file}.epoch;
    % choose config based on epoch
    config = get_behavior_config(session,epoch,config_path);
    
    % pulls up video frame and grabs coords
    coords_table = grab_coords(videoObj,session.behavioralTracking{1,file}.notes,config);
    
    % save data
    save_file = fullfile(basepath,[extractBefore(session.behavioralTracking{1,file}.notes,vid_type),'_maze_coords.csv']);
    writetable(coords_table,save_file);
    
end

end

function prompt = create_prompt(config)
%creates prompt so users can know which coordinate to collect

% use variables indicated inc config table
varnames = config.Properties.VariableNames;
prompt = []; % initialize prompt

% loop through and concatenate vaarname with column values
for i = 1:length(varnames)
    prompt = [prompt repmat(varnames(i),size(config,1),1) config.(varnames{i})];
end

% join to create cell array used in main
prompt = join(prompt);

end

function config_file = grab_coords(videoObj,vidname,config_path)
% Uses config to prompt users to define xy coordinates of an image
% (videoObj). Outputs xy coordinates in form of a table. 

% go to folder containing video & image
fig = figure;
imshow(readFrame(videoObj)) % display the first frame
axis('equal')
set(gcf,'CurrentCharacter','@')
hold on
i=1;

config_file = readtable(config_path);
prompt = create_prompt(config_file);

% let the user click around the coordinates
while true
    xlabel({prompt{i};'Press U key to redo coordinate'})
    [X,Y] = ginput(1);
    if isempty(X)
        break
    end
    %check if mistake made
    k = get(gcf,'CurrentCharacter');
    if k == 'u' % if user presses undo
        i = i - 1; % reset to previous iteration
        set(gcf,'CurrentCharacter','@')
        continue
    end
    coords(i,:)=[X,Y];
    scatter(coords(:,1),coords(:,2),'*r')
    i=i+1;
    if i > length(prompt)
        break
    end
end

% build csv
name = repmat(vidname,length(coords),1);
config_file.vid_name = name;
config_file.x = coords(:,1);
config_file.y = coords(:,2);


end

function config_name = get_behavior_config(session,epoch,config_path)
% uses tags in session.epoch.name to pull appropriate config saved in
% config_path. Currently limited to multi, object, context, and open_field.
% LB 2022
configs = dir(config_path);
configs = configs(~[configs.isdir]);
configs = {configs.name};

% look through session epochs
name = session.epochs{1, epoch}.name;

if contains(name,{'multi'})
    config_name = configs{find(contains(configs,'multi'))};
elseif contains(name,{'object'})
    config_name = configs{find(contains(configs,'object'))};
elseif contains(name,{'context','open_field'})
    config_name = configs{find(contains(configs,'open_field'))};
end

config_name = fullfile(config_path,config_name);

end