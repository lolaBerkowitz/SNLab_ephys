function df = get_ied_df(basepath)
% Load ied events for a given basepath and returns a dataframe consisting
% ofcolumns basepath, ied_peaks, ied_flag. 
%   basepath: basepath where IED data originates
%   ied_peaks: timestamp of when IED occurs
%   ied_flag: 1 = flagged as false detection, 0 = not flagged, likely real
%           IED. 
%
% Dataframe helpful for computing the accuracy of IED detection, but
% compiling IEDS and indicating if they were a false detection (flagged). 

basename = basenameFromBasepath(basepath);

% find the IED events 
events = dir(fullfile(basepath,'*.events.mat'));
ied_events = events(contains({events.name},{'interictal_spikes','IED'}));

if length(ied_events) == 1
	ied = load(fullfile(basepath,ied_events.name)) ;
else
    error('Multiple IED events found in basepath.')
end

fields = fieldnames(ied); 


 
ied_peaks = ied.(fields{1}).peaks;
if isfield(ied.(fields{1}),'flagged')
    ied_flag = ismember(ied.(fields{1}).peaks,ied.(fields{1}).peaks(ied.(fields{1}).flagged));
else
    ied_flag = repmat(NaN,length(ied_peaks),1); 
end
df = table;
df.basepath = repmat({basepath},length(ied_peaks),1); 
df.ied_peaks = ied_peaks;
df.ied_flag = ied_flag;

end