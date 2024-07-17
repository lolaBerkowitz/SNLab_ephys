function merged_events = merge_CE_events(event1,event2)
% Merge events takes two CellExplorer event files and merges them. For
% simplicity, flagged events are removed and added events are added from
% the individual event files before merging. 

% remove flagged events 
event1 = remove_flags(event1);
event2 = remove_flags(event2);

% include added events 
event1 = include_added_events(event1);
event2 = include_added_events(event2);

% consolidate the intervals from each event file
[intervals, ~] = ConsolidateIntervals([event2.timestamps;event1.timestamps]);


% gather peaks and keep ones that are included in the intervals
peaks = [event2.peaks(:); event1.peaks(:)];

[samples, ~] = Restrict(peaks,intervals,'sep',true);
merged_peaks = sort(samples{1});

% verify new peaks are in invervals 
if ~all(intervals(1) <= merged_peaks <= intervals(2))
    warning('Peaks are not within interval bounds. Suggest recalculating peaks.') 
elseif all(intervals(1) <= merged_peaks <= intervals(2))
    disp('Peaks validated')
end

% create new structure 
merged_events = struct; 
merged_events.timestamps = intervals; 
merged_events.peaks = sort(merged_peaks);
merged_events.duration = intervals(:,2) - intervals(:,1); 
merged_events.center = median(intervals,2);
merged_events.detectorinfo.detectorname = 'merge_CE_events';  
merged_events.detectorinfo.detectordate = date;

% rename struct to have the varname


end 



% Local functions below 
% remove flagged intervals 
function event_struct = remove_flags(event_struct)
if any(contains(fieldnames(event_struct),'flagged'))   
    event_struct.timestamps(event_struct.flagged,:) = [];
    event_struct.peaks(event_struct.flagged) = [];
end
end

% add intervals 
function event_struct = include_added_events(event_struct)
if any(contains(fieldnames(event_struct),'added'))
    event_struct.peaks = [event_struct.peaks(:); event_struct.added(:)];
end
end