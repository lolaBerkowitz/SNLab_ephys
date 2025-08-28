function optoStim = make_opto_stimulation_file(basepath,varargin)
%Make opto event files associated with stimulation and detection file loads 
% the processed digitalin.events file and retrieves the TTL pulses from
% the online detector. 

% input: 
% basepath: path (char/string) to data location of digitalin.events.mat
% file. 
% oupput: 
%  - optoStim and optoDetect mat files. These are CellExplorer formatted
%  events. 
%
p = inputParser;
addParameter(p,'stim_channel',6,@isnumeric) %snlab rig wiring for opto digitalin events
addParameter(p,'detect_channel',7,@isnumeric) %snlab rig wiring for opto digitalin events
addParameter(p,'detector','IED_60-80_Franksboard',@ischar) %snlab rig wiring for opto digitalin events

parse(p,varargin{:});

stim_channel = p.Results.stim_channel;
detection_channel = p.Results.detect_channel;
detector = p.Results.detector;

% load event file with stim 
load(fullfile(basepath,'digitalIn.events.mat'))

[~,basename] = fileparts(basepath);

% make optoStim
optoStim = struct;
optoStim.timestamps = [digitalIn.timestampsOn{stim_channel}, digitalIn.timestampsOff{stim_channel}];
optoStim.peaks = median(optoStim.timestamps,2);
optoStim.amplitude = zeros(length(optoStim.peaks),1);
optoStim.amplitudeUnits = 'AU';
optoStim.eventID = 1:length(optoStim.peaks);
optoStim.eventIDlabels = repmat({'LED_ON'},length(optoStim.peaks),1);
optoStim.eventIDbinary = false;
optoStim.center = median(optoStim.timestamps);
optoStim.duration = digitalIn.dur{stim_channel};
optoStim.detectorinfo = detector;

% make optoDetect
optoDetect = struct;
optoDetect.timestamps = [digitalIn.timestampsOn{detection_channel}, digitalIn.timestampsOff{detection_channel}];
optoDetect.peaks = median(optoDetect.timestamps,2);
optoDetect.amplitude = zeros(length(optoDetect.peaks),1);
optoDetect.amplitudeUnits = 'AU';
optoDetect.eventID = 1:length(optoDetect.peaks);
optoDetect.eventIDlabels = repmat({'Online_Detect'},length(optoDetect.peaks),1);
optoDetect.eventIDbinary = false;
optoDetect.center = median(optoDetect.timestamps);
optoDetect.duration = digitalIn.dur{stim_channel};
optoDetect.detectorinfo = detector;

% initialize the structure file 
save([basepath,filesep,basename,'.optoStim.events.mat'],'optoStim');
save([basepath,filesep,basename,'.optoDetect.events.mat'],'optoDetect');

end

