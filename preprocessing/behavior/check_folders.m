function check = check_folders(basepath)
% see if DLC files are present for each video

% find all .avi videos and dlc files in basepath
video_files = dir(fullfile(basepath,'*.avi'));
video_files = video_files(~cellfun(@(x) ismember(x(1,2),'._'), {video_files.name}));
video_files = extractBefore({video_files.name},'.avi');

dlc_files = get_dlc_files_in_basepath(basepath);
dlc_files = extractBefore({dlc_files.name},'DLC');

% see if each video has a DLC file
if sum(contains(video_files,dlc_files)) == length(video_files)
    check = true;
else
    check = false;
end


end