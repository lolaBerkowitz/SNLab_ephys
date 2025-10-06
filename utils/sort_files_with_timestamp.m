function file_struct = sort_files_with_timestamp(file_struct)
% Sorts files in dir output that contain a timestamp appended at the end of
% the file name (i.e. 3768N_pairing1_A-07122023090924.avi with sort by
% taking the digits after the - or 07122023090924. Returns sorted file
% structure. 

% LB 7/2023

% Extract the timestamps from the filenames
% timestamps = regexp({file_struct.name}, '-\d+', 'match');
% timestamps = cellfun(@(x) extractAfter(x,'-'),timestamps,'UniformOutput',false);
% timestamps = cellfun(@(x) str2double(x),timestamps,'UniformOutput',false);

names = {file_struct.name};

timestamps = [];
for i = 1:length(file_struct)
    name_i = names{i};
    
    if contains(name_i,'-0000')
        timestamp_i = regexp(name_i, '-\d+\-', 'match');
        timestamp_i = str2double(extractBetween(timestamp_i,'-','-'));
    else
        timestamp_i = regexp(name_i, '-\d+', 'match');
        timestamp_i = str2double(extractAfter(timestamp_i,'-'));
    end
    
    timestamps(i,1) = timestamp_i;
end

% Sort the file list based on the timestamps
[~, sortedIndices] = sort(timestamps);

% Sort the file list using the sorted indices
file_struct = file_struct(sortedIndices);
end