function make_events(varargin)

% make_Mergepointscreates a CellExplorer event file for each folder in a
% basepath, or for each video file in a basepath (behavior only studies). 

% input parser 
p = inputParser;
addParameter(p,'basepath',pwd,@isfolder); % by default, current folder
addParameter(p,'event_type','digitalIn',@isstring); % by default, current folder


parse(p,varargin{:});

basepath = p.Results.basepath; 
event_type = p.Results.event_type;
% See if subfolders exist

files = dir(basepath);
files = files(~ismember({files.name},{'.','..'}));

% Look through subfolders or identify videos to use 
if any([files.isdir])
    subfolders = files([files.isdir]);
end

% Obtain DLC file details 
dlc_files = find_dlc_files(basepath)

% Obtain 
switch event_type 
    case 'digitalIn'
        
    case 'MergePoints'
end
% Make the events.mat file that saves all the merge information
eventsfilename = fullfile(basepath,[basename,'.MergePoints.events.mat']);

MergePoints.timestamps = transitiontimes_sec; % in seconds
MergePoints.timestamps_samples = transitiontimes_samp; % in samples 
MergePoints.firstlasttimpoints_samples = firstlasttimepoints; 
MergePoints.foldernames = recordingnames; % names of videos 
MergePoints.filesmerged = datpaths; % basepaths of files merged 
MergePoints.filesizes = NaN; 
MergePoints.sizecheck = NaN; % 
MergePoints.detectorinfo.detectorname = 'bz_ConcatenateDats';
MergePoints.detectorinfo.detectiondate = datestr(now,'yyyy-mm-dd');

%Saving SleepStates
save(eventsfilename,'MergePoints');

eventsfilename = fullfile(basepath,[basename,'.MergePoints.events.mat']);

digitalIn.timestampsOn =  digital_on;
digitalIn.timestampsOff = digital_off;
digitalIn.ints =  digital_on;
digitalIn.dur = digital_on;
digitalIn.intsPeriods =  digital_on;

save(eventsfilename,'MergePoints');


end

function get_video_details(video_files)

for i = 1:length(video_files) 
% load video 
vidObj = VideoReader(fullfile(video_files.folder, video_files.name));


end

% find video files 
function dlc_files = find_dlc_files(basepath)

files = dir(basepath);
dlc_files = files(contains({files.name},'_filtered.csv'));

end

