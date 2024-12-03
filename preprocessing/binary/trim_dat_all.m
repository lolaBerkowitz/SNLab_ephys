% trim_dat_all allows the user to remove intervals from the binary files 
% (amplifier, auxiliary, digitalin, time) from Intan recordings.
% Useful in case of unplugs, bad wires, or other mishaps that happen during
% recordings. 

% Pipeline: 1) Copy origial files 2) Set intervals to remove 3) save new
% files with intervals removed. 

% Should be performed before any preprocessing. 

% Code adapted from ryanharvey1 cut_out_bad_intervals_from_session.m, edits
% by LB 2024


%% User defines intervals of time to remove (in samples interval_i by [start, end] ) 
bad_intervals = [],...
    [11612.2, 11890.2];

basepath = pwd;
basename = basenameFromBasepath(basepath);

%% load session for sample rate and make copy
load(fullfile(basepath,[basename, '.session.mat']))
fs = session.extracellular.sr;
save([basename, '.session_original.mat'], 'session')

%% Make copies of dat files 
movefile('auxiliary.dat', 'auxiliary_original.dat');
movefile('digitalin.dat', 'digitalin_original.dat');
% movefile('analogin.dat','analogin_original.dat');
movefile('time.dat', 'time_original.dat');
movefile(['amplifier', '.dat'], ['amplifier', '_original.dat']);

%% load memmap of amplifier and determin number of samples
d = dir(['amplifier', '_original.dat']);
nSamp = d.bytes / 2 / session.extracellular.nChannels;
data = memmapfile(['amplifier', '_original.dat'], 'Format', {'int16', [session.extracellular.nChannels, nSamp], 'x'});

session_domain = [0, (nSamp / fs) + 1];

%% locate good intervals
idx = ~InIntervals((1:nSamp) / fs, bad_intervals);
n_samples_new = sum(idx);

% get batch for amplifier (to reduce memory load)
batch = ceil(linspace(0,nSamp,10000));

%% Write files 

% use local function to write amplifier in a loop. 
write_amp(data,basepath,batch,idx)

% write auxiliary file
d = dir('auxiliary_original.dat');
nChannels = double(int64(d.bytes/2/nSamp));
data = memmapfile('auxiliary_original.dat', 'Format', {'uint16', [nChannels, nSamp], 'x'});
f = fopen('auxiliary.dat', 'w'); 
fwrite(f, data.Data.x(:, idx), 'uint16');
fclose(f);

% write digitalin file
d = dir('digitalin_original.dat');
nChannels = double(int64(d.bytes/2/nSamp));
data = memmapfile('digitalin_original.dat', 'Format', {'uint16', [nChannels, nSamp], 'x'});
f = fopen('digitalin.dat', 'w'); 
fwrite(f, data.Data.x(:, idx), 'uint16');
fclose(f);

% write time file
error('This part of the script does not reset time')
d = dir('time_original.dat');
nChannels = double(int64(d.bytes/2/nSamp));
data = memmapfile('time_original.dat', 'Format', {'int32', [nChannels, nSamp], 'x'});
f = fopen('time.dat', 'w'); 
fwrite(f, data.Data.x(:, idx), 'int32');
fclose(f);

%%%%%% Local function below

function write_amp(amp,basepath,batch,idx)
kept_samples = find(idx); 
amp_file = fopen(fullfile(basepath,'amplifier.dat'),'w');
% loop though batches
for ii = 1:length(batch)-1
    disp(['batch ',num2str(batch(ii)+1),' to ',num2str(batch(ii+1)),...
        '   ',num2str(ii),' of ',num2str(length(batch)-1)])
    % write to disk
    loop_idx = kept_samples(kept_samples >= batch(ii)+1 & kept_samples <= batch(ii+1));
    fwrite(amp_file,amp.Data.x(:,loop_idx), 'int16');
end
end