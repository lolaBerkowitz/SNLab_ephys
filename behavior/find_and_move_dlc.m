function find_and_move_dlc(dlc_folder,data_folder)
% find_and_move_dlc finds dlc output from dlc folder and moves to basepath
% located in data_folder.
%
% input:
% dlc_folder; folder where dlc output and videos are kept
% data_folder: path to subject folder or dataframe containing basepaths to data.
% output:
%


% Find all files from dlc_folder;
files = dir(dlc_folder);
files([files.isdir],:) = []; % remove directories

% load datafolder (handle input as csv of basepaths or directory)
if contains(data_folder,'.csv')
    df = readtable(data_folder);
else
    df = compile_sessions(data_folder);
end

basepaths = [df.basepath];
videos = files(contains({files.name},'.avi'));
videos = {videos.name}';
% Use data_folder to match where dlc files belong
for i = 1:length(basepaths)
    basepath = basepaths{i}{1};
    vid_name = extractBefore(videos,'.');
    
    % if video exists in basepath, transfer all files that contain video name
    if exist(fullfile(basepath,videos))
        
        temp = files(contains({files.name},vid_name));
        for ii = 1: length({temp.name})
            movefile(fullfile(temp(ii).folder,temp(ii).file),fullfile(basepath,temp(ii).file))
        end
    else
        continue
    end
    
    
end

% transfer dlc files
end