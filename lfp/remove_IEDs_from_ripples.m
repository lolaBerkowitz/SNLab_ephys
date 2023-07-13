function remove_IEDs_from_ripples(basepath)

basename = basenameFromBasepath(basepath);
% Load ripples, if ripple event overlaps with IED then define as
% IED ( Gelinas et al 2016)
load(fullfile(basepath,[basename,'.ripples.events.mat']))
load(fullfile(basepath,[basename,'.interictal_spikes_.events.mat']))

% Find IEDs that were also detected as physiological ripples
[status,interval,~] = InIntervals(interictal_spikes_.peaks,ripples.timestamps);

% Set those phyviological ripples as flagged
ripples.flagged = interval(status);

save(fullfile(basepath,[basename '.ripples.events.mat']),'ripples')
end