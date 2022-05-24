function get_MUA(varargin)
% get_MUA returns a CellExplorer timeseries container for events related to
% multiunit acitivity. 
%
% highpass filter from 300. Timestamps are
% selected from the z-scored envelope (anything above 3*SD of signal). 

% TO-DO: add median filter. As it is now, EMG artificats make up most of
% MUA. 

% LBerkowitz 2022

p = inputParser;
addParameter(p,'basepath',pwd,@isfolder) 
addParameter(p,'region','CA1',@ischar) 
addParameter(p,'channels',[],@isvec) 
addParameter(p,'high',300,@isdouble) 
addParameter(p,'sd_thresh',3,@isdouble)
addParameter(p,'binsz',.001,@isnumeric)
addParameter(p,'tSmooth',.016,@isnumeric);

parse(p,varargin{:});

basepath = p.Results.basepath;
region = p.Results.region;
high = p.Results.high;
sd_thresh = p.Results.sd_thresh;
binsz = p.Results.binsz;
tSmooth = p.Results.tSmooth;
channels = p.Results.channels;

% load session get recoring fs, 
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename,'.session.mat']))

fs = session.extracellular.sr;
n_channels = session.extracellular.nChannels;
bad_chan = session.channelTags.Bad.channels;

% finds channels used to compute MUA
if isfield(session,'brainRegions') & isempty(channels)
    fields = fieldnames(session.brainRegions);
    idx = contains(fieldnames(session.brainRegions),region);
    fields = fields(idx);
    for i = 1:length(fields)
        channels = [channels, session.brainRegions.(fields{i}).channels];
    end
    channels(ismember(channels,bad_chan)) = []; 
elseif ~isempty(channels)
    channels(ismember(channels,bad_chan)) = []; % remove bad channels
else
    disp('no brainRegions found or channels input found. Supply channels to compute MUA')
end

% find dat, obtain chunking info based on dat size 
contFile = fullfile(basepath,[basename,'.dat']);
file = dir(contFile);
samples = file.bytes/(n_channels * 2); %int16 = 2 bytes
amp = memmapfile(contFile,'Format',{'int16' [n_channels, samples] 'mapped'});
ts = linspace(0,samples/fs,samples); % estimate time stamps in seconds 

% initalize high pass filter 
[b1, a1] = butter(3, high/(fs/2), 'high');

% Filter each channel
tic
fprintf('finding MUA on channel...');
mua_spk = [];
for ii = channels
    fprintf(' %d ', ii);
    
    % high pass filter data with forwards backwards filter and get envelope
    datr = abs(filtfilt(b1, a1, double(amp.Data.mapped(ii,:))));
    
    % Find times where voltage surpases SD of sigal 
    mua_idx = datr > sd_thresh*std(datr);
    
    % convert index to ts and save 
    mua_spk = [mua_spk, ts(mua_idx)];
end

% bin over time
bin_ts = 0:binsz:(samples/fs);
spkhist = histcounts(mua_spk, bin_ts);

% smooth binned spikes
fsize = tSmooth/binsz;
gfilt = fspecial('gauss', [10*fsize 1], fsize);
spkhist = conv(spkhist, gfilt, 'same');

% save as timeseries container
timeseries = struct;
timeseries.data = spkhist;
timeseries.timestamps = bin_ts;
timeseries.precision = 'double';
timeseries.units = 'spike_counts';
timeseries.nchannels = 1;
timeseries.channelNames = 1;
timeseries.sr = 1/binsz;
timeseries.nsamples = length(spkhist);
timeseries.description = ['binned MUA for area',region];
timeseries.processinginfo.function = 'get_MUA ';
timeseries.processinginfo.date = datetime;
timeseries.processinginfo.sourceFileName = [basename,'.dat'];
toc
save(fullfile(basepath,[basename,'.',region,'_MUA.timeseries.mat']),'timeseries')

% s = 100000;
% figure;
% subplot(2,1,1)
% plot(abs(datr(1:s))); hold on; plot(1:s,repmat(60,s,1),'r');
% subplot(2,1,2)
% plot(amp.Data.mapped(ii,1:s),'b');

% zscore and find threshold positive and negative crossings

% extract timestamps of crossing 

% compile results 




end