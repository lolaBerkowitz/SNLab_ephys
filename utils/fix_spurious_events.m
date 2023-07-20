% Manual approach to fix events coming from unused digitalin inputs.
% Sometimes the relays discharge leading to spurious events on all event
% channels. This manually removes events on used from unused channels
basepath = pwd;
unused_channels = [2,5,8:13];
event_channels = [3,4];

load(fullfile(basepath,'digitalIn.events.mat'))
%% get unique events from unused channels 
spurious_events = unique(cat(1,digitalIn.timestampsOn{1,unused_channels}));
for channel = event_channels
    idx = ismember(digitalIn.timestampsOn{1,channel},spurious_events);
    digitalIn.timestampsOn{1,channel}(idx) = [];
    digitalIn.timestampsOff{1,channel}(idx) = [];
    disp(['Do these time differences make sense? for event channel : ',num2str(channel)])
    diff(digitalIn.timestampsOn{1,channel})/60
end

%% Check that the events make sense realtive to notes


%% Remove data from unused channels 
for unused = unused_channels
    digitalIn.timestampsOn{1,unused} = [];
    digitalIn.timestampsOff{1,unused} = [];
end

%% set off relative to on 
for channel = event_channels
    digitalIn.timestampsOff{1,channel} = digitalIn.timestampsOn{1,channel}+.1;
end

%% Save digitalIn back to basepath 

save([basepath,filesep,'digitalIn.events.mat'],'digitalIn');


