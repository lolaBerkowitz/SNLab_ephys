function digitalIn = make_events(varargin)

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

% Obtain
switch event_type
    case 'digitalIn'
        
        % initialize
        eventsfilename = fullfile([basepath,filesep,'digitalin.events.mat']);
        % create
        digitalIn = make_digitalIn_from_dlc(files);
        % save
        save(eventsfilename,'digitalIn');
    case 'MergePoints'
        disp('TO-DO as of 7/13/2023')
        %         % Make the events.mat file that saves all the merge information
        %         eventsfilename = fullfile(basepath,[basename,'.MergePoints.events.mat']);
        %
        %         MergePoints.timestamps = transitiontimes_sec; % in seconds
        %         MergePoints.timestamps_samples = transitiontimes_samp; % in samples
        %         MergePoints.firstlasttimpoints_samples = firstlasttimepoints;
        %         MergePoints.foldernames = recordingnames; % names of videos
        %         MergePoints.filesmerged = datpaths; % basepaths of files merged
        %         MergePoints.filesizes = NaN;
        %         MergePoints.sizecheck = NaN; %
        %         MergePoints.detectorinfo.detectorname = 'bz_ConcatenateDats';
        %         MergePoints.detectorinfo.detectiondate = datestr(now,'yyyy-mm-dd');
        %
        %         %Saving SleepStates
        %         save(eventsfilename,'MergePoints');
end



end

function [fs,duration] =  get_video_details(vid_file_path)

% load video
vid_obj = VideoReader(vid_file_path);
fs = vid_obj.FrameRate;
duration = vid_obj.Duration;

end


% make digitalIn events from DLC
function digitalIn = make_digitalIn_from_dlc(files)

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

% sort
dlc_files = sort_files(dlc_files);

% loop through files and generate digitalIn
if length(dlc_files) > 1
    % loop through files
    ts = [];
    offset = [];
    for i = 1:length(dlc_files)
        dlc_file = fullfile(dlc_files(i).folder,dlc_files(i).name);
        if isempty(ts)
            temp_ts = get_ts_from_dlc(dlc_file,files); %timestamps
        else
            temp_ts = get_ts_from_dlc(dlc_file,files) + max(ts); %add offset
        end
        ts = [ts temp_ts];
        digitalIn.timestampsOn{1,2}(i,1) = min(temp_ts);
        digitalIn.timestampsOff{1,2}(i,1) = max(temp_ts);
    end
    digitalIn.timestampsOn{1,1} = ts';
    digitalIn.timestampsOff{1,1} = ts';
    
else isempty(dlc_files)
    digitalIn.timestampsOn{1,1} = nan;
    digitalIn.timestampsOff{1,1} = nan;
    digitalIn.timestampsOn{1,2} = nan;
    digitalIn.timestampsOff{1,2} = nan;
end


digitalIn = make_digitalin_vars(digitalIn);

end

function file_struct = sort_files(file_struct)
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

function ts = get_ts_from_dlc(dlc_file,files)
dlc_file_name = basenameFromBasepath(dlc_file);
vid_name = extractBefore(dlc_file_name,'DLC_');
all_vid_files = files(contains({files.name},{'.avi','.mpg'}));
all_vid_files(cellfun(@(x) ismember(x(1,2),'._'), {all_vid_files.name})) = [];
current_vid = all_vid_files(contains({all_vid_files.name},vid_name));
vid_file_path = fullfile(current_vid.folder,current_vid.name);
%%
[fs,~] =  get_video_details(vid_file_path);
%%
df = load_dlc_csv(dlc_file);
col_names = fieldnames(df);
samples = 1:length(df.(col_names{1})); % get length from first column of df
ts = linspace(1/fs,length(samples)/fs,length(samples));
end

function digitalIn = make_digitalin_vars(digitalIn)

%% and recalculate other variables
for ii = 1:size(digitalIn.timestampsOn,2)
    lag = 100;
    if isempty(digitalIn.timestampsOn{ii})
        continue
    end
    % intervals
    d = zeros(2,max([size(digitalIn.timestampsOn{ii},1) size(digitalIn.timestampsOff{ii},1)]));
    d(1,1:size(digitalIn.timestampsOn{ii},1)) = digitalIn.timestampsOn{ii};
    d(2,1:size(digitalIn.timestampsOff{ii},1)) = digitalIn.timestampsOff{ii};
    
    if isempty(d)
        continue
    end
    
    if d(1,1) > d(2,1)
        d = flip(d,1);
    end
    if d(2,end) == 0; d(2,end) = nan; end
    digitalIn.ints{ii} = d;
    digitalIn.dur{ii} = digitalIn.ints{ii}(2,:) - digitalIn.ints{ii}(1,:); % durantion
    
    clear intsPeriods
    intsPeriods(1,1) = d(1,1); % find stimulation intervals
    intPeaks =find(diff(d(1,:))>lag);
    for jj = 1:length(intPeaks)
        intsPeriods(jj,2) = d(2,intPeaks(jj));
        intsPeriods(jj+1,1) = d(1,intPeaks(jj)+1);
    end
    intsPeriods(end,2) = d(2,end);
    digitalIn.intsPeriods{ii} = intsPeriods;
end
end
