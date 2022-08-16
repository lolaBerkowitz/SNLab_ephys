function get_OL_XY(varargin)
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
%  add functionality to add data to CellExplorer mat files
%  add ability to scroll through frames for start/end trials
%
% LB 2022

% input parser
p = inputParser;
p.addParameter('basepath',pwd,@isfolder);
p.addParameter('vid_type','.avi',@ischar);
p.addParameter('vid_time',300,@isnumeric); % time of video to load in seconds
p.addParameter('ephys',true,@islogical); % updates basename.animal.behavior file with coords
p.addParameter('multi_chamber',false,@islogical); % if there are multiple open fields/tasks in same video
p.addParameter('config_path','C:\Users\schafferlab\github\SNLab_ephys\behavior\behavior_configs',@islogical); % time of video to load in seconds


p.parse(varargin{:})
basepath = p.Results.basepath; % not used currently
vid_type = p.Results.vid_type; % not used currently
vid_time = p.Results.vid_time; % not used currently
ephys = p.Results.ephys; % not used currently
multi_chamber = p.Results.multi_chamber;
config_path = p.Results.config_path;

% find all video files indicated by vid_type
vid_files = dir(fullfile(basepath,['*',vid_type]));

if isempty(vid_files)
    error('No videos found. Check video type and path')
end

% get basename for loading CellExplorer containers 
basename = basenameFromBasepath(basepath);

if ephys % saves csv and appends coordinates to general behavior file
    % load config file 
    config = readtable(fullfile(config_path,'ephys_config.csv'));
    
    % grab coordinates
    main(basepath,vid_files,config,vid_type,vid_time)
       
    % check if behavior file exsists and if not make one
    % load session file
    
    % check if file was already made and if so load
    if exist([basepath,filesep,[basename,'.animal.behavior.mat']],'file')
        load([basepath,filesep,[basename,'.animal.behavior.mat']],'behavior');
    else % generate one
        general_behavior_file_SNlab(basepath)
    end
    
    % load session file
    session = loadSession(basepath,basename);
    
    % get coords 
    main(basepath,vid_files,config,vid_type,vid_time)

    % load coord csv 
%     coords = dir(fullfile(basepath,'*object_maze_coords.csv'));
    
elseif ~ephys % time of video to load in seconds
    % load config
    config = readtable(fullfile(config_path,'ephys_config.csv'));
    
    main(basepath,vid_files,config,vid_type,vid_time)
    
elseif ~ephys && multi_chamber
    
    config = readtable(fullfile(config_path,'multi_chamber_config.csv'));
    
    main(basepath,vid_files,config,vid_type,vid_time)
end

end

function main(basepath,vid_files,config,vid_type,vid_time)
% runs main process of looping through videos in basepath, pulling up
% image via local grab_coords function which allows user to collect coordinate data
% and outputs coords for object A center, object A edge, object B center,
% object B edge, corner center A-D (inner corner of maze). 
% 

% LB 2022

basename = basenameFromBasepath(basepath);

% loop through video
for file = 1:length(vid_files) %loop through folders containing subject videos
    
    vid_path = fullfile(basepath,vid_files(file).name);
    
    %Pull up video
    videoObj = VideoReader(vid_path,'CurrentTime',vid_time); % load video starting at vid_time
    
    % pulls up video frame and grabs coords
    coords_table = grab_coords(videoObj,vid_files(file).name,config);
    
    % save data
    save_file = fullfile(basepath,[basename,'_',extractBefore(vid_files(file).name,vid_type),'object_maze_coords.csv']);
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

function config = grab_coords(videoObj,vidname,config)

% go to folder containing video & image
imshow(readFrame(videoObj)) % display the first frame
set(gcf,'CurrentCharacter','@')
hold on
i=1;

prompt = create_prompt(config);
% let the user click around the coordinates
while true
    title(prompt{i})
    xlabel('Press U key to redo coordinate')
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
config.name = repmat(vidname,length(coords),1);
config.x = coords(:,1);
config.y = coords(:,2);

end