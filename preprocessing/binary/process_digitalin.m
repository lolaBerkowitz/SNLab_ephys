function digitalIn = process_digitalin(basepath,fs,varargin)
% code adapted from getDigitalin.m in neurocode
% (https://github.com/ayalab1/neurocode/tree/master/preprocessing)

p = inputParser;
addParameter(p,'dat_name','digitalin.dat',@ischar) %snlab rig wiring for events
addParameter(p,'filter_channels',[1 2 3 4 5 6 7 8 9 10 11 12 13 14 15],@isnumeric) %channels to apply median filter in case of noise
addParameter(p,'lag',20,@isnumeric)
parse(p,varargin{:});

dat_name = p.Results.dat_name;
filter_channels = p.Results.filter_channels;
lag = p.Results.lag; 

% % check if digitalin.events.map exists and load if so
if isfile([basepath,filesep,'digitalin.events.mat'])
    load([basepath,filesep,'digitalin.events.mat'])
    return 
end
    

contFile = fullfile(basepath,dat_name);
D.Data = memmapfile(contFile,'Format','uint16','writable',false);

% intan has 16 digitalin channels, set n channel and n channel + 1 (Nchan2)
% to account for zero and 1 based indexing. 
Nchan = 16;
Nchan2 = 17;
for k = 1:Nchan
    test = bitand(D.Data.Data,2^(Nchan-k)) > 0;
    
    % filter by 5 to remove single sample noise events
    if any((filter_channels-1)==(Nchan-k))
        test = medfilt1(single(test),5) > 0; 
    end

    digital_on{Nchan2-k} = find(diff(test) == 1);
    digital_off{Nchan2-k} = find(diff(test) == -1);
end


% take pulses and create output
for ii = 1:size(digital_on,2)
    if ~isempty(digital_on{ii})
        % take timestamp in seconds
        digitalIn.timestampsOn{ii} = digital_on{ii}/fs;
        digitalIn.timestampsOff{ii} = digital_off{ii}/fs;
        
        % intervals
        d = zeros(2,max([size(digitalIn.timestampsOn{ii},1) size(digitalIn.timestampsOff{ii},1)]));
        d(1,1:size(digitalIn.timestampsOn{ii},1)) = digitalIn.timestampsOn{ii};
        d(2,1:size(digitalIn.timestampsOff{ii},1)) = digitalIn.timestampsOff{ii};
        if d(1,1) > d(2,1)
            d = flip(d,1);
        end
        if d(2,end) == 0; d(2,end) = nan; end
        digitalIn.ints{ii} = d;
        %duration of intervals
        digitalIn.dur{ii} = digitalIn.ints{ii}(2,:) - digitalIn.ints{ii}(1,:); % 
        
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
    
end
if ~exist('digitalIn','var')
    disp('no events detected, saving empty struct')
    digitalIn.timestampsOn =  digital_on;
    digitalIn.timestampsOff = digital_off;
    digitalIn.ints =  digital_on;
    digitalIn.dur = digital_on;
    digitalIn.intsPeriods =  digital_on;
end
save([basepath,filesep,'digitalIn.events.mat'],'digitalIn');
end