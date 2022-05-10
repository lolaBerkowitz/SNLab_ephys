function digitalIn = process_digitalin(basepath,fs,varargin)
% code adapted from getDigitalin.m in neurocode
% (https://github.com/ayalab1/neurocode/tree/master/preprocessing)

p = inputParser;
addParameter(p,'dat_name','digitalin.dat',@ischar) %snlab rig wiring for events
parse(p,varargin{:});

dat_name = p.Results.dat_name;

% check if digitalin.events.map exists and load if so
if isfile([basepath,filesep,'digitalin.events.mat'])
    load([basepath,filesep,'digitalin.events.mat'])
    return 
end
    
lag = 20; % This pertains to a period for known event (in this case stimulation period).
%           keeping for now but will need to update if we ever have events
%           for stimulation. Making large cause events are based on double-clicks
%           of varying length LB 2/22


contFile = fullfile(basepath,dat_name);
% file = dir(fullfile(basepath,dat_name));
% samples = file.bytes/2/16; % 16 is n_channels for intan digitial in
D.Data = memmapfile(contFile,'Format','uint16','writable',false);

digital_word2 = double(D.Data.Data);
Nchan = 16;
Nchan2 = 17;
for k = 1:Nchan
    tester(:,Nchan2-k) = (digital_word2 - 2^(Nchan-k))>=0;
    digital_word2 = digital_word2 - tester(:,Nchan2-k)*2^(Nchan-k);
    test = tester(:,Nchan2-k) == 1;
    test2 = diff(test);
    pulses{Nchan2-k} = find(test2 == 1);
    pulses2{Nchan2-k} = find(test2 == -1);
    data(k,:) = test;
end
digital_on = pulses;
digital_off = pulses2;


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
    
end
if ~exist('digitalIn','var')
    disp('no events detected, saving empty struct')
    digitalIn.timestampsOn =  digital_on;
    digitalIn.timestampsOff = digital_off;
end
save([basepath,filesep,'digitalIn.events.mat'],'digitalIn');
end