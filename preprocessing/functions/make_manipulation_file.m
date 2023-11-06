function manipulation = make_manipulation_file(basepath,varargin)
%Make manipulation file loads the analogin.dat file and retirieves the
%channel associated with optogenetic stimulation waveform. 
p = inputParser;
addParameter(p,'stim_channel',16,@isnumeric) %snlab rig wiring for opto digitalin events
addParameter(p,'event_label','LED_ON',@ischar) %snlab rig wiring for events

parse(p,varargin{:});

stim_channel = p.Results.stim_channel;
event_label = p.Results.event_label;
% load event file with stim 
load(fullfile(basepath,'digitalIn.events.mat'))

basename = basenameFromBasepath(basepath);
manipulation = struct;
manipulation.timestamps = [digitalIn.timestampsOn{stim_channel}, digitalIn.timestampsOff{stim_channel}];
manipulation.peaks = median(manipulation.timestamps);
manipulation.amplitude = zeros(length(manipulation.peaks),1);
manipulation.amplitudeUnits = 'AU';
manipulation.eventID = 1:length(peaks);
manipulation.eventIDlabels = repmat({event_label},length(manipulation.peaks),1);
manipulation.eventIDbinary = false;
manipulation.center = median(manipulation.timestamps);
manipulation.duration = digitalIn.dur{stim_channel};
manipulation.detectorinfo = 'SWRDetector v.00';

% initialize the structure file 
save([basepath,filesep,basename,'.opto_stimulation.manipulation.mat'],'manipulation');
end

