basepath = pwd
dat_files = dir(fullfile(basepath,'**','*amplifier.dat'));

recordingnames = {dat_files.folder}
names2sort = cellfun(@(X) str2num(X(end-5:end)),recordingnames,'UniformOutput',false);
names2sort = cell2mat(names2sort);
[~,I] = sort(names2sort);
recordingnames = recordingnames(I);

dat_files = dat_files(I);

basename = basenameFromBasepath(basepath);

sessionInfo = LoadXml(fullfile(basepath,[basename, '.xml']));
nSamp = [];
for didx = 1:length(dat_files)
    % use dir to get size of dat file in bytes
    % determine number of bytes per sample
    dataTypeNBytes = numel(typecast(cast(0, 'int16'), 'uint8'));
    % Number of samples per channel
    nSamp(didx) = dat_files(didx).bytes/(sessionInfo.nChannels*dataTypeNBytes);
end
cumsum_nSamp = cumsum([nSamp]);
starts = [0,cumsum_nSamp(1:end-1)];
transitiontimes_samp = [starts',cumsum_nSamp'];
transitiontimes_sec = transitiontimes_samp./sessionInfo.SampleRate;

firstlasttimepoints = [zeros(length(nSamp),1),nSamp'];

recordingnames = [];
for didx = 1:length(dat_files)
    [filepath,name,ext] = fileparts(dat_files(didx).folder);
    recordingnames{1,didx} = name;
end

MergePoints.timestamps = transitiontimes_sec;
MergePoints.timestamps_samples = transitiontimes_samp;
MergePoints.firstlasttimpoints_samples = firstlasttimepoints;
MergePoints.foldernames = recordingnames;
MergePoints.detectorinfo.detectorname = 'custom';
MergePoints.detectorinfo.detectiondate = datestr(now,'yyyy-mm-dd');

save(fullfile(basepath,[basename,'.MergePoints.events.mat']),'MergePoints');