
function preprocess_session(basepath,varargin)
% Function for processing electrophysiology data in the SNLab adapted from
% preprocessSession.m ayalab/neurocode repo
%
% Dependencies
%   - external packages: neurocode,KS2Wrapper
%   - schafferlab/github/Kilosort
%   - schafferlab/github/CellExplorer
%   - schafferlab/github/neurocode/utils
%
% Data organization assumptions
%   - Data should be organized with each recordding session saved within an animal folder i.e. 'animal_id/basename'
%   - 'basename' is used by CellExplorer and other buzcode functions.
%     Therefore, it is important to make this name unique.
%   - XML with channel mapping completed should be done and located in the
%     basepath.

% LBerkowitz 2021

p = inputParser;
addParameter(p,'analogInputs',false,@islogical);
addParameter(p,'multishank',false,@islogical);
addParameter(p,'analogChannels',[],@isnumeric);
addParameter(p,'digitalInputs',false,@islogical);
addParameter(p,'digitalChannels',[],@isnumeric);
addParameter(p,'getAcceleration',true,@islogical);
addParameter(p,'concatenateDat',false,@islogical);
addParameter(p,'getLFP',true,@islogical);
addParameter(p,'getEMG',true,@islogical);
addParameter(p,'stateScore',true,@islogical);
addParameter(p,'tracking',true,@islogical);
addParameter(p,'kilosort',false,@islogical);
addParameter(p,'removeNoise',false,@islogical);
addParameter(p,'ssd_path','C:\kilo_temp',@ischar);    % Path to SSD disk.
addParameter(p,'lfp_fs',1250,@isnumeric);
addParameter(p,'specialChannels',[30,59,17],@isnumeric) % default for ASYpoly2A64
addParameter(p,'acquisition_event_flag',false,@islogical) % if you have start and stop recording
% events else default uses first and last event as indicies for start and stop of recording
addParameter(p,'check_epochs',true,@islogical) % fixes bouncy events
addParameter(p,'maze_size',30,@isnumeric); % in cm


parse(p,varargin{:});

analogInputs = p.Results.analogInputs;
multishank = p.Results.multishank;
digitalInputs = p.Results.digitalInputs;
getAcceleration = p.Results.getAcceleration;
getLFP = p.Results.getLFP;
getEMG = p.Results.getEMG;
concatenateDat = p.Results.concatenateDat;
kilosort = p.Results.kilosort;

% cleanArtifacts = p.Results.cleanArtifacts;
stateScore = p.Results.stateScore;
tracking = p.Results.tracking;
ssd_path = p.Results.ssd_path;
lfp_fs = p.Results.lfp_fs;
specialChannels = p.Results.specialChannels;
check_epochs = p.Results.check_epochs;

% Prepare dat files and prepare metadata

% Set basename from folder
[~,basename] = fileparts(basepath);

% If non-contiguous recordings were taken, merge the dat files into one
% session. Default is false.
if concatenateDat
    concatenate_dats(basepath);
end


% rename amplifier.dat to basename.dat
if ~isempty(dir([basepath,filesep,'amplifier.dat']))
    disp(['renaming amplifer.dat to ',basename,'.dat'])
    % create command
    command = ['rename "',basepath,filesep,'amplifier.dat"',' ',basename,'.dat'];
    system(command); % run through system command prompt
end

% lets also rename the xml if present.
if ~isempty(dir([basepath,filesep,'amplifier.xml']))
    disp(['renaming amplifer.xml to ',basename,'.xml'])
    % create command
    command = ['rename "',basepath,filesep,'amplifier.xml"',' ',basename,'.xml'];
    system(command); % run through system command prompt
end

% Create SessionInfo (only run this once!)
session = sessionTemplate(basepath,'showGUI',false);
save([basepath,filesep,[basename,'.session.mat']],'session')


% Load rhd file 
% Load info.rhd for recording parameters
[~, ~, ~, ~,~, ~, frequency_parameters,~] = read_Intan_RHD2000_file_snlab(basepath);

% Process additional inputs
% Analog inputs
% check the two different fucntions for delaing with analog inputs and proably rename them
if analogInputs
    % analog pulses ...
    [pulses] = getAnalogPulses('samplingRate',frequency_parameters.board_adc_sample_rate);
end

% Digital inputs
if digitalInputs
    
    digitalIn = process_digitalin(basepath,frequency_parameters.amplifier_sample_rate);

end

% Epochs derived from digital inputs for multianimal recordings
update_epochs('basepath',basepath,'annotate',check_epochs)

% Auxilary input
if getAcceleration
    accel = computeIntanAccel('saveMat',true);
end

% get lfp
if getLFP
    % create downsampled lfp low-pass filtered lfp file
    LFPfromDat(basepath,'outFs',lfp_fs);
end

% Calcuate estimated emg
if getEMG
    if isempty(specialChannels)
        EMGFromLFP = getEMGFromLFP(basepath,'overwrite',false,'noPrompts',true,...
            'rejectChannels',session.channelTags.Bad.channels);
    else
        
        EMGFromLFP = getEMGFromLFP(basepath,'overwrite',false,'noPrompts',true,...
            'specialChannels',specialChannels,'rejectChannels',session.channelTags.Bad.channels);
    end
end

% Get brain states
% an automatic way of flaging bad channels is needed
if stateScore
    
    SleepScoreMaster(basepath,'rejectChannels',session.channelTags.Bad.channels,'overwrite',false); % takes lfp in base 0
    
end


if kilosort

    % For Kilosort: create channelmap
    if ~isempty(session.channelTags.Bad.channels)
        createChannelMapFile_KSW(basepath,basename,'staggered',session.channelTags.Bad.channels);
    else
         createChannelMapFile_KSW(basepath,basename,'staggered');
    end

    % creating a folder on the ssd for chanmap,dat, and xml
    ssd_folder = fullfile(ssd_path, basename);
    mkdir(ssd_folder);
    
    % Copy chanmap,basename.dat, and xml
    disp('Copying basename.dat, basename.xml, and channelmap to ssd')
    
    disp('Saving dat file to ssd')
    command = ['robocopy "',basepath,'" ',ssd_folder,' ',basename,'.dat'];
    system(command);
    
    disp('Saving xml to ssd')
    command = ['robocopy "',basepath,'" ',ssd_folder,' ',basename,'.xml'];
    system(command);
    
    disp('Saving channel_map to ssd')
    command = ['robocopy "',basepath,'" ',ssd_folder,' chanMap.mat'];
    system(command);
    
    % Spike sort using kilosort 1 (data on ssd)
    run_ks1(basepath,'ssd_folder',ssd_folder)
end
% Get tracking positions

if tracking
    % update epochs so process_tracking can run
    gui_session
    % saves general behavior file, restrictxy, maze_coords, and updates
    % session.behavioralTracking. Requires dlc output in basepath
    process_tracking(basepath)
    
end
    
end


