function file_struct = sort_files_with_timestamp(file_struct)
% Sorts files in dir output that contain a timestamp appended at the end of
% the file name (i.e. 3768N_pairing1_A-07122023090924.avi with sort by
% taking the digits after the - or 07122023090924. Returns sorted file
% structure. 

% LB 7/2023

% Extract the timestamps from the filenames
timestamps = regexp({file_struct.name}, '-\d+', 'match');
timestamps = cellfun(@(x) extractAfter(x,'-'),timestamps,'UniformOutput',false);
timestamps = cellfun(@(x) str2double(x),timestamps,'UniformOutput',false);

% Convert the timestamps to numeric values

% Sort the file list based on the timestamps
[~, sortedIndices] = sort( [timestamps{:}]);

% Sort the file list using the sorted indices
file_struct = file_struct(sortedIndices);
end