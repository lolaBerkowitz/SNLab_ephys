function optoStim = make_optoStim_file(basepath,varargin)
%Make manipulation file loads the analogin.dat file and retirieves the
%channel associated with optogenetic stimulation waveform. 
p = inputParser;
addParameter(p,'stim_channel',7,@isnumeric) %snlab rig wiring for opto digitalin events
addParameter(p,'detector','Franksboard_v4',@ischar) %snlab rig wiring for opto digitalin events
addParameter(p,'event_label','LED_ON',@ischar) %snlab rig wiring for events

parse(p,varargin{:});

stim_channel = p.Results.stim_channel;
event_label = p.Results.event_label;
detector = p.Results.detector;

% load event file with stim 
load(fullfile(basepath,'digitalIn.events.mat'))

[~,basename] = fileparts(basepath);
optoStim = struct;
optoStim.timestamps = [digitalIn.timestampsOn{stim_channel}, digitalIn.timestampsOff{stim_channel}];
optoStim.peaks = median(optoStim.timestamps,2);
optoStim.amplitude = zeros(length(optoStim.peaks),1);
optoStim.amplitudeUnits = 'AU';
optoStim.eventID = 1:length(optoStim.peaks);
optoStim.eventIDlabels = repmat({event_label},length(optoStim.peaks),1);
optoStim.eventIDbinary = false;
optoStim.center = median(optoStim.timestamps);
optoStim.duration = digitalIn.dur{stim_channel};
optoStim.detectorinfo = detector;


% initialize the structure file 
save([basepath,filesep,basename,'.optoStim.events.mat'],'optoStim');
end

