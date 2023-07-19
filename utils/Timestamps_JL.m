
% Get all the videos from the directory
basepath='/Volumes/sn data server 3/laura_berkowitz/behavior_validation/appps1_cpp/data/cohort2/3779/3779_pretest_pairing_day01';
file_struct = dir(fullfile(basepath,'*.avi'));
fil
% maybe use later
File_name=basepath(1+max(strfind(basepath,'/')):end); 
% Extract the timestamps from the filenames
timestamps = regexp({file_struct.name}, '-\d+', 'match');
timestamps = cellfun(@(x) extractAfter(x,'-'),timestamps,'UniformOutput',false);
timestamps = cellfun(@(x) str2double(x),timestamps,'UniformOutput',false);

% Convert the timestamps to numeric values

% Sort the file list based on the timestamps
[~, sortedIndices] = sort( [timestamps{:}]);

% convert to matrix so you can apply matrix functions 
timestamps = cell2mat(timestamps); 
%convert to just HHMMSS
timestamps=mod(timestamps,1000000);
hours=floor(timestamps/10000);
minutes=floor((timestamps-hours*10000)/100);
seconds=mod(timestamps,100);
%new matrix of time of day in seconds
newtimestamps=hours*3600+60*minutes+seconds;
%Sort and find the difference
newtimestamps=sort(newtimestamps,'descend');
difference=abs(diff(newtimestamps/60));
row={File_name};
for i=1:1:length(file_struct)-1
    row=[row,difference(i)]
end
% load video
vid_obj = VideoReader('/Volumes/sn data server 3/laura_berkowitz/behavior_validation/appps1_cpp/data/cohort2/3769/3769_pretest_pairing_day01/3769L_pairing1_A-07122023111233.avi');
fs = vid_obj.FrameRate;
duration = vid_obj.Duration;






% Given the format of the timestamps is MMDDYYYYHHMMSS, extract the time
% information so you can compute the inter-video interval or time between
% the start time of each video. 