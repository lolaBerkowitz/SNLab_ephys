function dlc_files = get_dlc_files_in_basepath(basepath)

% use dir to find all files in basepath
files = dir(basepath);
files = files(~ismember({files.name},{'.','..'})); % remove non-files

% get csv files
csv_files = files(contains({files.name},'.csv'));

% and keep only those with DLC in name
dlc_files = csv_files(contains({csv_files.name},'DLC'));

% if dlc file was filtered, keep that one, else keep unfiltered
toss_idx = [];
for i = 1:length(dlc_files)
    cur_dlc = extractBefore(dlc_files(i).name,'DLC');
    if sum(contains({dlc_files.name},cur_dlc)) > 1
        % keep filtered
        toss_idx = [toss_idx; find(contains({dlc_files.name},cur_dlc) & ~contains({dlc_files.name},'filtered.csv'))];
    end
end
% remove duplicate csv
dlc_files(toss_idx,:) = [];
end