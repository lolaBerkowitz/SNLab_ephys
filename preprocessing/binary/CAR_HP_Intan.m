function CAR_HP_Intan(dat_path, varargin)
% Loads intan amplifier.dat and performs CAR to prepare for spike sorting. 
%
% Subtracts median of each channel, then subtracts median of each time
% point and also high-pass filters (default at 300 Hz)
% Does so in chunks, users buffers to avoid artefacts at edges
%
% filename should be the complete path to an intan.rhd file
% The same directory should contain the 'amplifier.dat' file that is the
% raw data


% from
% https://github.com/PaulMAnderson/IntanTools/blob/master/CAR_HP_Intan.m
% adapted heavily for SNLab_ephys by LBerkowitz 08/2021

%% unpack inputs 
p = inputParser;
p.addParameter("runFilter",true)
p.addParameter("high_pass",300)
p.parse(varargin{:});

runFilter = p.Results.runFilter;
hi_pass = p.Results.high_pass;


% get header file info and load
d = dir([dat_path,filesep,'*.rhd']);

intanRec = read_Intan_RHD2000_file_snlab([d.folder,filesep],d.name);


% get amplifier info
if ~exist([d.folder,filesep,'amplifier.dat'],'file')
    error('amplifier.dat file not in header directory...')
else  
    amp_info = dir([d.folder,filesep,'amplifier.dat']);
end

%% Setup Parameters
% should make chunk size as big as possible so that the medians of the
% channels differ little from chunk to chunk.
chunkSize  = 2^20;
bufferSize = 2^10;

numChannels = length(intanRec.amplifier_channels);
numSamples  = amp_info.bytes/numChannels/2; % samples = bytes/channels/2 (2 bits per int16 sample)
fs       = intanRec.frequency_parameters.amplifier_sample_rate;
numChunks   = ceil(numSamples./chunkSize);


filename = fullfile(amp_info.folder,amp_info.name);
if runFilter
    outFilename = fullfile(amp_info.folder,[amp_info.name(1:end-4) '_CAR_HP.dat']); % file name: amplifier_CAR_HP.dat
else
    outFilename = fullfile(amp_info.folder, [amp_info.name(1:end-4) '_CAR.dat']);
end

fid    = fopen(filename,'r');
fidOut = fopen(outFilename,'w');

amplifierMap = memmapfile(filename,...
'Format', {
'int16', [numChannels numSamples], 'data'
});

% initalize high pass filter 
if runFilter
    % Create high pass filter
    [b1, a1] = butter(3, hi_pass/fs*2, 'high');
end
% filteredData = zeros(numChannels,numSamples,'int16');

%% Loop through the chunks
for chunkI = 1:numChunks
    tic     
    fprintf('Loading Chunk %d of %d... \n',chunkI,numChunks)
    if chunkI == 1
        startPoint = 1;
        endPoint   = chunkSize;
        chunk = amplifierMap.Data.data(:,1:chunkSize+bufferSize);
        chunk = [zeros(numChannels,bufferSize,'int16') chunk];
    
    elseif chunkI == numChunks
        startPoint = (chunkSize * (chunkI-1)) + 1;
        endPoint   = numSamples;
        chunk      =  amplifierMap.Data.data(:,...
            chunkSize * (chunkI-1) + 1 - bufferSize : numSamples);
        lastChunkSize = size(chunk,2);
        if lastChunkSize < chunkSize + 2 * bufferSize
            chunk = [chunk zeros(numChannels, ...
            (chunkSize + 2 * bufferSize) - lastChunkSize,'int16')];
        end
    else
        chunk = amplifierMap.Data.data(:,...
            chunkSize * (chunkI-1) + 1 - bufferSize : ...
             chunkSize*chunkI  + bufferSize);
        startPoint = (chunkSize * (chunkI-1)) + 1;
        endPoint   = chunkSize * (chunkI);
    end
        
    
    %% Baseline, Rereference and Filter
    fprintf('Re-Referencing...\n')
    chunk = bsxfun(@minus, chunk, median(chunk,2)); % subtract median of each channel
    chunk = bsxfun(@minus, chunk, median(chunk,1)); % subtract median of each time point
      
    if runFilter
        fprintf('Filtering...\n')
        % high pass filter data with forwards backwards filter
        chunk = filtfilt(b1, a1, double(chunk));
%         chunk = eegfilt(chunk,fs, loFreq, hiFreq);
    end
    %% Write out the data
    
    if chunkI == numChunks 
        fwrite(fidOut,chunk(:,bufferSize+1:lastChunkSize),'int16');   
    else    
        fwrite(fidOut,chunk(:,bufferSize+1:end-bufferSize),'int16');    
    end
    % Optionally combine the data here

       % filteredData(:,startPoint:endPoint) = ...
       %    chunk(:, bufferSize+1:end-bufferSize);

end
toc
fclose(fid);
fclose(fidOut);