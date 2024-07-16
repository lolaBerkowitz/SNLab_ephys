function manipulation = make_manipulation_file(basepath,varargin)
%Make manipulation file loads the analogin.dat file and retirieves the
%channel associated with optogenetic stimulation waveform. 
p = inputParser;
addParameter(p,'stim_channel',16,@isnumeric) %snlab rig wiring for opto digitalin events
addParameter(p,'ramp_channel',15,@isnumeric) %snlab rig wiring for opto digitalin events

addParameter(p,'event_label','LED_ON',@ischar) %snlab rig wiring for events

parse(p,varargin{:});

stim_channel = p.Results.stim_channel;
event_label = p.Results.event_label;
ramp_channel = p.Results.ramp_channel;
% load event file with stim 
load(fullfile(basepath,'digitalIn.events.mat'))

basename = basenameFromBasepath(basepath);
manipulation = struct;
manipulation.timestamps = [digitalIn.timestampsOn{ramp_channel}, digitalIn.timestampsOff{ramp_channel}];
manipulation.peaks = median(manipulation.timestamps,2);
manipulation.amplitude = zeros(length(manipulation.peaks),1);
manipulation.amplitudeUnits = 'AU';
manipulation.eventID = 1:length(manipulation.peaks);
manipulation.eventIDlabels = repmat({event_label},length(manipulation.peaks),1);
manipulation.eventIDbinary = false;
manipulation.center = median(manipulation.timestamps);
manipulation.duration = digitalIn.dur{ramp_channel};
manipulation.detectorinfo = 'SWRDetector v4';

% initialize the structure file 
save([basepath,filesep,basename,'.opto_stimulation.events.mat'],'opto_stim');
end

