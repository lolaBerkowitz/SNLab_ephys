function debounce_events(basepath,digitalIn_ch)
% debound_events finds timestamps that are at least 1 second apart for
% digitalIn channel and saves back to digitalIn.events.mat
% 
% load event file
load(fullfile(basepath,['digitalIn.events.mat']))

try
events_on = digitalIn.timestampsOn{1, digitalIn_ch}; 
catch
    disp('channel index not valid.') 
    return
end

events_off = digitalIn.timestampsOff{1, digitalIn_ch}; 

% Find gaps in time that are larger than 1 second 
idx_off = find(diff(events_on)  > 1);
idx_on = idx_off - 1; %for east just take timestamp before

% Save last timestamp before the jump, and append final timestamp in the
% end
digitalIn.timestampsOn{1,digitalIn_ch} = digitalIn.timestampsOn{1, digitalIn_ch}([idx_on; length(events_on)-1]);
digitalIn.timestampsOff{1,digitalIn_ch} = digitalIn.timestampsOff{1, digitalIn_ch}([idx_off; length(events_on)]);
digitalIn.ints = [digitalIn.timestampsOn{1,digitalIn_ch} digitalIn.timestampsOff{1,digitalIn_ch}];
digitalIn.dur = [digitalIn.timestampsOn{1,digitalIn_ch} - digitalIn.timestampsOff{1,digitalIn_ch}];

save([basepath,filesep,'digitalIn.events.mat'],'digitalIn');

end