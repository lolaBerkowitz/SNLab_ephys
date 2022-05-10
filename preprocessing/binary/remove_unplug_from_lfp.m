basepath = 'D:\app_ps1\data\hpc04\hpc04_day09_220221_094101';
basename = basenameFromBasepath(basepath);
epoch = [7268.332 7800.357]; % epoch to remove in seconds

% locate lfp file
fileName = fullfile(basepath,[basename,'.lfp']);

% use dir to get size of dat file in bytes
filenamestruct = dir(fileName);
dataTypeNBytes = numel(typecast(cast(0, 'int16'), 'uint8')); % determine number of bytes per sample
nSamp = filenamestruct.bytes/(64*dataTypeNBytes);  % Number of samples per channel
mmf = memmapfile(fileName, 'Format', {'int16', [64,nSamp], 'Data'},'Writable',true);

%load lfp to get timestamps
lfp = getLFP('all','basepath',basepath);

% find index to remove
idx = lfp.timestamps >= epoch(1) & lfp.timestamps <= epoch(2);
mmf.Data.Data(:,idx) = 0;

clear mmf lfp





