function opto_stimulation = make_opto_stimulation_file(basepath,varargin)
%Make manipulation file loads the analogin.dat file and retirieves the
%channel associated with optogenetic stimulation waveform. 
p = inputParser;
addParameter(p,'stim_channel',16,@isnumeric) %snlab rig wiring for opto digitalin events
addParameter(p,'detector','IED_60-80_Franksboard',@ischar) %snlab rig wiring for opto digitalin events
addParameter(p,'event_label','LED_ON',@ischar) %snlab rig wiring for events

parse(p,varargin{:});

stim_channel = p.Results.stim_channel;
event_label = p.Results.event_label;
detector = p.Results.detector;

% load event file with stim 
load(fullfile(basepath,'digitalIn.events.mat'))

[~,basename] = fileparts(basepath);
opto_stimulation = struct;
opto_stimulation.timestamps = [digitalIn.timestampsOn{stim_channel}, digitalIn.timestampsOff{stim_channel}];
opto_stimulation.peaks = median(opto_stimulation.timestamps,2);
opto_stimulation.amplitude = zeros(length(opto_stimulation.peaks),1);
opto_stimulation.amplitudeUnits = 'AU';
opto_stimulation.eventID = 1:length(opto_stimulation.peaks);
opto_stimulation.eventIDlabels = repmat({event_label},length(opto_stimulation.peaks),1);
opto_stimulation.eventIDbinary = false;
opto_stimulation.center = median(opto_stimulation.timestamps);
opto_stimulation.duration = digitalIn.dur{stim_channel};
opto_stimulation.detectorinfo = detector;

% initialize the structure file 
save([basepath,filesep,basename,'.opto_stimulation.events.mat'],'opto_stimulation');
end

