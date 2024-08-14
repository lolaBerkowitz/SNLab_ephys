% check_events processes digitalin events from intan digitalin.dat. Allows
% user to verify number of events and duration of epochs in case of event
% errors that may arise during acquisition. 
%
% Load info.rhd for recording parameters
[amplifier_channels, ~, aux_input_channels, ~,...
    ~, ~, frequency_parameters,board_adc_channels ] = ...
    read_Intan_RHD2000_file_snlab(basepath);

% load exsisting events
digitalIn = process_digitalin(basepath,frequency_parameters.amplifier_sample_rate);

%% remove pulses from electrical noise that is present on all channels 
unused_channel = [5,8,9,10,11,12,13];
used_channel = [2,3,4];

for channel = used_channel
    for unused = unused_channel
        digitalIn.timestampsOff{1,channel}(ismember(digitalIn.timestampsOff{1,channel},digitalIn.timestampsOff{1,unused})) = [];
        digitalIn.timestampsOn{1,channel}(ismember(digitalIn.timestampsOn{1,channel},digitalIn.timestampsOn{1,unused})) = [];
    end
end

% remove unused ttl pulses
for unused = unused_channel
    digitalIn.timestampsOff{1,unused} = [];
    digitalIn.timestampsOn{1,unused} = [];
end

% Examine digitialIn structure and verify events match written notes
% if edited, delete ints, dur, and intsPeriods and remake belo
%%

digitalIn.timestampsOff{1,2} = [];
digitalIn.timestampsOff{1,3} = [];
digitalIn.timestampsOff{1,4} = [];


digitalIn.timestampsOff{1,2} = digitalIn.timestampsOn{1,2} + .3;
digitalIn.timestampsOff{1,3} = digitalIn.timestampsOn{1,3} + .3;
digitalIn.timestampsOff{1,4} = digitalIn.timestampsOn{1,4} + .3;




%% and recalculate other variables
for ii = 1:size(digitalIn.timestampsOn,2)
    lag = 100;
    if isempty(digitalIn.timestampsOn{ii})
        continue
    end
    % intervals
    d = zeros(2,max([size(digitalIn.timestampsOn{ii},1) size(digitalIn.timestampsOff{ii},1)]));
    d(1,1:size(digitalIn.timestampsOn{ii},1)) = digitalIn.timestampsOn{ii};
    d(2,1:size(digitalIn.timestampsOff{ii},1)) = digitalIn.timestampsOff{ii};
    
    if isempty(d)
        continue
    end
    
    if d(1,1) > d(2,1)
        d = flip(d,1);
    end
    if d(2,end) == 0; d(2,end) = nan; end
    digitalIn.ints{ii} = d;
    digitalIn.dur{ii} = digitalIn.ints{ii}(2,:) - digitalIn.ints{ii}(1,:); % durantion
    
    clear intsPeriods
    intsPeriods(1,1) = d(1,1); % find stimulation intervals
    intPeaks =find(diff(d(1,:))>lag);
    for jj = 1:length(intPeaks)
        intsPeriods(jj,2) = d(2,intPeaks(jj));
        intsPeriods(jj+1,1) = d(1,intPeaks(jj)+1);
    end
    intsPeriods(end,2) = d(2,end);
    digitalIn.intsPeriods{ii} = intsPeriods;
end
save([basepath,filesep,'digitalIn.events.mat'],'digitalIn');


