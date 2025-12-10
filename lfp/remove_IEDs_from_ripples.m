function remove_IEDs_from_ripples(basepath)
% Load ripples, if ripple event overlaps with IED then define as
% IED ( Gelinas et al 2016)
% LB 2023

basename = basenameFromBasepath(basepath);

rip_file = fullfile(basepath,[basename,'.ripples.events.mat']);
ied_file = fullfile(basepath,[basename,'.interictal_spikes_offline.events.mat']);

% Load ripples and ieds
if isfile(rip_file) & isfile(ied_file)
    load(fullfile(basepath,[basename,'.ripples.events.mat']),'ripples')
    load(fullfile(basepath,[basename,'.interictal_spikes_offline.events.mat']),'interictal_spikes_offline')
else
    disp('Missing Ripples or IEDs, skipping')
    return
end

% remove flagged events (sometimes these are ripples
if isfield(interictal_spikes_offline,'flagged')
    interictal_spikes_offline.peaks(interictal_spikes_offline.flagged) = [];
end

% set interval bounds - anything 10ms before IED peak and 200ms after 
bounds = [interictal_spikes_offline.peaks' - 0.010, interictal_spikes_offline.peaks' + 0.2];

% Find IEDs that were also detected as physiological ripples
[status,~,~] = InIntervals(ripples.peaks,bounds);

% Set those phyviological ripples as flagged
if isfield(ripples,'flagged')
    ripples.flagged = unique(sort([ripples.flagged; find(status)]));
else
    ripples.flagged = find(status);
end

save(fullfile(basepath,[basename '.ripples.events.mat']),'ripples')
end