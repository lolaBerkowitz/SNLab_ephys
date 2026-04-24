[f, p] = uigetfile('amplifier.dat');
cd(p)

fs = 30e3;

nChannels = 64;
threshold = .3;

min_time_motionless = 10; % seconds
time = 0;
iter = 0;
rd_win = 1*fs;
res_fs = fs/rd_win;
curr_time = 0;
pass = logical([]);
acc_values = [];
downsample_factor = 16;
maxfreqband = floor(max([625, fs / downsample_factor / 2]));
xcorr_freqband = [275, 300, maxfreqband - 25, maxfreqband]; % Hz

filename = fullfile(p, f);

fprintf('Processing file: %s\n', filename);

xcorr_chs = pick_channels(filename);

while (1)
    % read aux signal
    data=[];
    while(size(data,1)<rd_win)
        data = readmulti_frank([p, f], nChannels, xcorr_chs, curr_time, curr_time+rd_win, 'int16');
    end
    data_ds = downsample(data, downsample_factor); % Downsample by factor

    % Filter the data
    xcorr_freqband = [275, 300, maxfreqband - 25, maxfreqband]; % Hz
    filtered_data = filtsig_in(data_ds, fs/downsample_factor, xcorr_freqband);
    
    % Calculate the EMG value
    emg_value = computeEMGValue(filtered_data');
    % disp(emg_value)

    curr_time = curr_time + size(data,1);
    % filter and get acceleration
    % acc_fil = filter(b, a, data);
    % acc = abs(acc_fil(:, end));
    % acc_value = mean(acc(:));
    % if below threshold, count as sleeping
    % disp(acc_value)
    if (emg_value < threshold)
        
        time = time + 1;
    end
    % acc_values = [acc_values;acc_value];
    pass = [pass;emg_value < threshold];
    % display every 30 seconds
    if iter > 30*res_fs
        % multiple by 120 because reading file is ~120 times
        %   slower than sample rate
        % disp_time((time)/res_fs)
        iter = 0;
        % iter_time = toc;
        % disp(['iter time: ',num2str(iter_time)])
        % tic
        pass_intervals = findIntervals(pass);
        durations = (pass_intervals(:,2) - pass_intervals(:,1)) / res_fs;
        disp_time(sum(durations(durations >= min_time_motionless)))
    end
    iter = iter + 1;
end

function out = findIntervals(booString)
%findIntervals - find start and end indices from boolstring
booString = reshape(booString,[1,length(booString)]);
starts = find(diff([false booString])>0);
ends = find(diff([booString false])<0);

out = [starts', ends'];
end

% pretty display
function disp_time(x_time)
% calc total duration
obj_duration = seconds(x_time);
if obj_duration < seconds(1)
    duration_str = datestr(obj_duration, 'FFF');
    units = 'ms';
elseif obj_duration < seconds(60)
    duration_str = datestr(obj_duration, 'SS.FFF');
    units = 'seconds';
elseif obj_duration < seconds(3600)
    duration_str = datestr(obj_duration, 'MM:SS.FFF');
    units = 'minutes';
elseif obj_duration < seconds(86400)
    duration_str = datestr(obj_duration, 'HH:MM:SS.FFF');
    units = 'hours';
else
    duration_str = datestr(obj_duration, 'DD:HH:MM:SS.FFF');
    units = 'days';
end
fprintf('sleep time %s %s \n', ...
    ...
    duration_str, ...
    units);
end


function processEMGFile(filename, nChannels, sample_rate,...
    downsample_factor, chunk_size, threshold, maxfreqband, output_file)
% PROCESS_EMG_FILE Processes an actively acquired file in chunks, filtering data and calculating EMG values.
%
% Parameters
% ----------
% filename : string
%     Path to the data file being actively acquired.
% nChannels : int
%     Number of channels in the data.
% sample_rate : double
%     Sampling rate of the amplifier in Hz.
% downsample_factor : int
%     Factor by which to downsample the data before filtering.
% filter_cutoff : double
%     Cutoff frequency for the band-pass filter in Hz.
% chunk_size : int
%     Number of samples to process per iteration (before downsampling).
%
% Returns
% -------
% None. Results are displayed in the console.

% Open the file and initialize loop variables
fprintf('Processing file: %s\n', filename);

xcorr_chs = pick_channels(filename);

while true
    try
        % Read a chunk of data using readmulti_frank
        raw_data = readmulti_frank(filename, nChannels, xcorr_chs, chunk_size, 0, 'int16');

        % Downsample the data
        data_ds = downsample(raw_data, downsample_factor); % Downsample by factor

        % Filter the data
        xcorr_freqband = [275, 300, maxfreqband - 25, maxfreqband]; % Hz
        filtered_data = filtsig_in(data_ds, sample_rate/downsample_factor, xcorr_freqband);

        % Calculate the EMG value
        emg_value = computeEMGValue(filtered_data');
        disp(emg_value)

        if (emg_value > threshold)
            value = 0;
        else
            value = 1;
        end
        try
            fh = fopen(output_file, 'w');
            fwrite(fh, value, 'int32');
            fclose(fh);
        catch
            warning("can't write file")
        end
    catch ME
        % Handle potential errors (e.g., file access issues or data reading issues)
        fprintf('Error processing chunk %s\n', ME.message);
        pause(1); % Pause briefly before retrying
    end
end
end

function xcorr_chs = pick_channels(filename)

%% Pick channels and to analyze
% get spike groups,
% pick every other one... unless specialshanks, in which case pick non-adjacent
%This is potentially dangerous in combination with rejectChannels... i.e.
%what if you pick every other shank but then the ones you pick are all
%reject because noisy shank.

basepath = fileparts(filename);
% xcorrs_chs is a list of channels that will be loaded
% spkgrpstouse is a list of spike groups to find channels from
sessionInfo = LoadXml(fullfile(basepath, 'amplifier.xml'));

channels = [sessionInfo.AnatGrps.Channels];
skipped_channel_idx = ([sessionInfo.AnatGrps.Skip]);
rejectChannels = channels(skipped_channel_idx == 1);

SpkGrps = {sessionInfo.SpkGrps.Channels};


% get list of spike groups (aka shanks) that should be used
usablechannels = [];
spkgrpstouse = [];
if length(SpkGrps) > 1
    n = 1;
else
    n = 5;
end
for gidx = 1:length(SpkGrps)
    usableshankchannels{gidx} = setdiff(SpkGrps{gidx}, rejectChannels);
    usablechannels = cat(2, usablechannels, usableshankchannels{gidx});
    if ~isempty(usableshankchannels{gidx})
        spkgrpstouse = cat(2, spkgrpstouse, gidx);
    end
end

% okay, lets limit usableshankchannels with spkgrpstouse. Thus, empty shanks will be excluded.
usableshankchannels = usableshankchannels(spkgrpstouse);

% get list of channels (1 from each good spike group)
xcorr_chs = [];
for gidx = 1:length(usableshankchannels)
    %grab random channel on each shank
    if ~isempty(usableshankchannels{gidx})
        randChfromShank = usableshankchannels{gidx}(randi(length(usableshankchannels{gidx}), n, 1));
        xcorr_chs = [xcorr_chs, randChfromShank];

    end
end
xcorr_chs = unique(xcorr_chs);
end

function emg_value = computeEMGValue(filtered_data)
% COMPUTE_EMG_VALUE Calculates the average correlation across channels.
%
% Parameters
% ----------
% filtered_data : matrix
%     Filtered data matrix of size [nChannels x nSamples].
%
% Returns
% -------
% emg_value : double
%     Average correlation across channels.

nChannels = size(filtered_data, 1);
corr_sum = 0;
count = 0;

% Iterate over all pairs of channels to compute pairwise correlations
for i = 1:nChannels
    for j = i + 1:nChannels
        % corr_sum = corr_sum + corrcoef(filtered_data(i, :)', filtered_data(j, :)');
        correlation = corrcoef(filtered_data(i, :)', filtered_data(j, :)');
        corr_sum = corr_sum + correlation(1,2);
        count = count + 1;
    end
end
emg_value = corr_sum / count;
end

function [filt_sig, Filt] = filtsig_in(sig, Fs, filtband_or_Filt)
% [filt_sig, Filt] = filtsig(sig, dt_ms, filtband_or_Filt)
%
% Created by: Erik Schomburg, 2011

if isnumeric(filtband_or_Filt)

    h = fdesign.bandpass( ...
        filtband_or_Filt(1), ...
        filtband_or_Filt(2), ...
        filtband_or_Filt(3), ...
        filtband_or_Filt(4), ...
        60, ...
        1, ...
        60, ...
        Fs);

    Filt = design(h, 'butter', 'MatchExactly', 'passband');
else
    Filt = filtband_or_Filt;
end

if ~isempty(sig)
    if iscell(sig)
        filt_sig = cell(size(sig));
        for i = 1:length(sig(:))
            filt_sig{i} = filter(Filt, sig{i});
            filt_sig{i} = filter(Filt, filt_sig{i}(end:-1:1));
            filt_sig{i} = filt_sig{i}(end:-1:1);
        end
    elseif ((size(sig, 1) > 1) && (size(sig, 2) > 1))
        filt_sig = zeros(size(sig));
        for i = 1:size(filt_sig, 2)
            filt_sig(:, i) = filter(Filt, sig(:, i));
            filt_sig(:, i) = filter(Filt, filt_sig(end:-1:1, i));
            filt_sig(:, i) = filt_sig(end:-1:1, i);
        end
    else
        filt_sig = filter(Filt, sig);
        filt_sig = filter(Filt, filt_sig(end:-1:1));
        filt_sig = filt_sig(end:-1:1);
    end
else
    filt_sig = [];
end
end

function [xml, rxml] = LoadXml(fbasename, varargin)
%function [xml, rxml] = LoadXml(FileBase)
%
% loads the xml file using xmltools (have to have it in the path)
% rxml returns it's original layout - very messy structure but contains all
% the xml file contents.
% xml - is the ouput structure which is backwards compatible to LoadPar
% output, so you can use it instead ..also loads some usefull stuff -
% Anatomoical groups with Skips , Spike electrode groups
% more can be added later (e.g. parameters of the process scripts)
% this script is written for xml version 1.1 .. older version doesn't work.
% additions are welcome
xml = struct;

xmli = strfind(fbasename, '.xml');
if isempty(xmli)
    fbasename = [fbasename, '.xml'];
end
rxml = xmltools(fbasename);
try
    rxml = rxml.child(2);
catch
    rxml = rxml.child;
end
% from this level all children are the different parameters fields
xml.FileName = fbasename;

for i = 1:length(rxml.child)

    switch lower(rxml.child(i).tag)

        case 'generalinfo'
            xml.Date = rxml.child(i).child(1).value; % date of xml file creation?

        case 'acquisitionsystem'
            xml.nBits = str2num(rxml.child(i).child(1).value); % number of bits of the file
            xml.nChannels = str2num(rxml.child(i).child(2).value);
            xml.SampleRate = str2num(rxml.child(i).child(3).value);
            xml.SampleTime = 1e6 / xml.SampleRate; %to make backwards compatible
            xml.VoltageRange = str2num(rxml.child(i).child(4).value);
            xml.Amplification = str2num(rxml.child(i).child(5).value);
            xml.Offset = str2num(rxml.child(i).child(6).value);

        case 'fieldpotentials'
            xml.lfpSampleRate = str2num(rxml.child(i).child.value);

        case 'anatomicaldescription'
            tmp = rxml.child(i).child.child;
            for grpI = 1:length(tmp)
                for chI = 1:length(tmp(grpI).child)
                    xml.AnatGrps(grpI).Channels(chI) = str2num(tmp(grpI).child(chI).value);
                    xml.AnatGrps(grpI).Skip(chI) = str2num(tmp(grpI).child(chI).attribs.value);
                end
            end

        case 'spikedetection'
            if ~isempty(rxml.child(i).child)
                tmp = rxml.child(i).child.child;
                for grpI = 1:length(tmp)
                    for chI = 1:length(tmp(grpI).child(1).child)
                        xml.SpkGrps(grpI).Channels(chI) = str2num(tmp(grpI).child(1).child(chI).value);
                    end
                    if length(tmp(grpI).child) > 1
                        xml.SpkGrps(grpI).nSamples = str2num(tmp(grpI).child(2).value);
                        xml.SpkGrps(grpI).PeakSample = str2num(tmp(grpI).child(3).value);
                        xml.SpkGrps(grpI).nFeatures = str2num(tmp(grpI).child(4).value);
                    end
                    %backwards compatibility
                    xml.nElecGps = length(tmp);
                    xml.ElecGp{grpI} = xml.SpkGrps(grpI).Channels;
                end
            else
                xml.nElecGps = 0;
            end


        case 'programs'
            tmp = rxml.child(i).child;
            for i = 1:length(tmp)
                if strcmp(tmp(i).child(1).value, 'process_mhipass')
                    for j = 1:length(tmp(i).child(2).child)
                        if strcmp(tmp(i).child(2).child(j).child(1).value, 'frequency')
                            xml.HiPassFreq = str2num(tmp(i).child(2).child(j).child(2).value);
                            break
                        end
                    end
                end
            end
    end
end

if ~isfield(xml, 'AnatGrps')
    xml.AnatGrps.Channels = [];
    xml.AnatGrps.Skip = [];
end
end

function [eeg, fb] = readmulti_frank(fname, numchannel, chselect, read_start, read_until, precision, b_skip)
% eeg is output, fb is size of file in bytes
% Reads multi-channel recording file to a matrix
% last argument is optional (if omitted, it will read all the
% channels.
%
% From the Buzsaki lab (obtained 4/5/2010).
% revised by zifang zhao @ 2014-5-1 increased 2 input to
% control the range of file reading
read_start = round(read_start);
read_until = round(read_until);
if nargin < 6 %precision and skip
    precision = 'int16';
end
if nargin < 7 %skip
    b_skip = 0;
end
fileinfo = dir(fname);
if nargin == 2
    datafile = fopen(fname, 'r');
    eeg = fread(datafile, [numchannel, inf], 'int16');
    fclose(datafile);
    eeg = eeg';
    return
end
fb = fileinfo(1).bytes;
byte_len = get_format_bytes(precision);
numel_all = floor(fb/byte_len/numchannel);
fb = numel_all * byte_len * numchannel;
if nargin >= 3
    % the real buffer will be buffersize * numch * 2 bytes
    % (int16 = 2bytes)
    if nargin < 4
        read_until = numel_all;
    end
    buffersize = 4096;
    % get file size, and calculate the number of samples per channel
    if read_start < 0
        read_start = read_start + numel_all - 1;
        if read_until == 0
            read_until = numel_all;
        end
    end

    read_start(read_start < 0) = 0;
    read_until = read_until + 1;
    read_until(read_until > numel_all) = numel_all;
    read_start_byte = read_start * byte_len * numchannel;
    read_until_byte = read_until * byte_len * numchannel;
    numel = read_until - read_start;


    eeg = zeros(numel, length(chselect));

    %% original method
    numel1 = 0;
    %  numelm=0;
    datafile = fopen(fname, 'r');
    state = fseek(datafile, read_start_byte, 'bof');
    %  while ~feof(datafile),
    while ftell(datafile) < read_until_byte && ~feof(datafile) && state == 0
        len_left = read_until_byte - ftell(datafile);
        if len_left >= buffersize * numchannel * byte_len
            [data, count] = fread(datafile, [numchannel, buffersize], precision, b_skip); %can be improved,vectorize,arrayfun,multi-threading, zifangzhao@4.24
        else
            [data, count] = fread(datafile, [numchannel, ceil(len_left/numchannel/byte_len)], precision, b_skip); %can be improved,vectorize,arrayfun,multi-threading, zifangzhao@4.24
        end
        numelm = (count / numchannel); %numelm = count/numchannel;
        if numelm > 0
            eeg(numel1+1:numel1+numelm, :) = data(chselect, :)';
            numel1 = numel1 + numelm;
        end
    end
end
fclose(datafile);
end

function len = get_format_bytes(format)
len = 2;
end