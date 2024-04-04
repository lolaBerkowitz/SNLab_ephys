
% Procedure for processing electrophysiology data in the SNLab
%
% Dependencies
%   - schafferlab/github/Kilosort
%   - schafferlab/github/CellExplorer
%   - schafferlab/github/neurocode
%
% Data organization assumptions
%   - Data should be organized with each recordding session saved within an animal folder i.e. 'animal_id/basename'
%   - 'basename' is used for general data organization. Therefore, it is important to make this name unique.
%
%
% Pipeline Overview: 
%
%  0.  Split-dat (if multiple animals are recorded in same sesion)
%  1.  Run DLC from anaconda prompt (env: DLC-GPU)
%  2.  Preprocess (create lfp, get digital signals, behavior tracking, kilosort)
%  3.  Manual Curation in Phy 
%  4.  Compile and save in CellExplorer data structure

% LBerkowitz 2021

%% paths for individual research projects 
% laura = 'Y:\laura_berkowitz\VR_ephys\data';
laura = 'Y:\laura_berkowitz\alz_stim\data';
% laura = 'Y:\laura_berkowitz\app_ps1_ephys\data';

%% Multi-animal recording session 
% For sessions that record from multiple headstages from separate animals.
% For SNLab, assumes one animal per active port (64 channel electrodes) as
% of 2/22
% subject folder corresponds to port A, B, C, D, respectively
subject_order = {'grigri','etrier',[],[]};

% folder where dat files reside that need to be split
data_path ='Y:\laura_berkowitz\alz_stim\data\to_split\grigri_day17_etrier_day20_240320_124022';

% project folder where subjects data should be saved
save_path = {laura,laura,[],[]}; 
    
% split dat files a
split_dat(data_path,save_path, subject_order,'trim_dat',false)

%% Process tracking (Done first so tracking can be updated in preprocess_session)
% Open anaconda prompt and open gui
% enter commands: 
% activate conda DLC-GPU
% ipython
% import deeplabcut
% deeplabcut.launch_dlc() 
% Load config for video to analyze (ie.
% 30cm_open_field-berkowitz-2022-07-14 for 30cm open fields)
% Once project is loaded, navigate to analyze_video tab and load the videos
% to be analyzed. Let them run overnight while split dat runs. 

%% Single Session Preprocess
basepath = 'Y:\laura_berkowitz\alz_stim\data\etrier\etrier_day17_231130_112605';

% Check xml/channel mapping: verify channel map, skip bad channels, and save 
% make_xml(basepath)

% Preprocess (create lfp, kilosort)


% OPTO ONLY - Preprocess (create lfp, kilosort)
preprocess_session(basepath,'digitalInputs',true,'kilosort',true,'tracking',false)


% multiple shank 
preprocess_session(basepath,'digitalInputs',false,'kilosort',true,'tracking',true,'specialChannels',[])


%% Batch preprocess (must make sure xml is made/accurate before running)
subject_folder = 'Y:\laura_berkowitz\app_ps1_ephys\data\hpc09'; %subject main folder (i.e. ~\data\hpc01)

preprocess_batch(subject_folder)

%% Clean up kilo results in Phy
disp ('Curate the kilosort results in Phy ')

% In anaconda prompt, cd to kilosort folder.
% cd ks2_folder i.e. cd C:\kilo_temp\2021-09-16_test_210916_135552
% conda activate phy2
% phy template-gui params.py

%% Copy back over to data file
disp ('Manually copy the kilosort folder from the ssd_path to the main data folder.')

%% Extract spike times and waveforms for sorted clusters
cd(basepath)
basename = basenameFromBasepath(basepath);
% update 
session = loadSession(basepath,basename);
gui_session

% run channel mapping to update session with 
channel_mapping('basepath',basepath)
gui_session

%% find ripples
% first input [ripple channel, sharp wave channel] using intan channels. 
ripple_channel = 62%session.brainRegions.CA1sp.channels(end-1);
sharp_wave_channel = 58%session.brainRegions.CA1sr.channels(end-3);
noise_channel = 46;
ripples = DetectSWR([ripple_channel,...
    sharp_wave_channel,...
    noise_channel],...
    'thresSDrip',[.25,.5],...
    'thresSDswD',[.25,.5],...
    'saveMat',true,...
    'check',true,...
    'forceDetect',true);

%% Compute basic cell metrics
cell_metrics = ProcessCellMetrics('session',session,'manualAdjustMonoSyn',false);
% GUI to manually curate cell classification. This is not neccesary at this point.
% It is more useful when you have multiple sessions
cell_metrics = CellExplorer('metrics',cell_metrics);
channel_mapping('basepath',basepath)


%% To load multiple sessions into CellExplorer
% After you have multiple sessions spike sorted and processed, you can open
% them all in cell explorer.

% assuming your current directory is (../project/animal/session), this will
% cd to project level (for example: A:\Data\GirardeauG)
cd('..\..')

% make this current directory into variable
data_path = pwd;

% look for all the cell_metrics.cellinfo.mat files
% files = dir([data_path,'\**\*.cell_metrics.cellinfo.mat']);
basepath = [];
basename = [];
% pull out basepaths and basenames
% for i = 1:length(files)
%     basepath{i} = files(i).folder;
%     basename{i} = basenameFromBasepath(files(i).folder);
% end

objectlocationpaths = readtable("Y:\laura_berkowitz\app_ps1_ephys\analysis\object_location_paths.csv");

for i = 1:length(objectlocationpaths.basepath)
    basepath{i} = objectlocationpaths.basepath{i};
    basename{i} = basenameFromBasepath(objectlocationpaths.basepath{i});
end

% load all cell metrics
cell_metrics = loadCellMetricsBatch('basepaths',basepath,'basenames',basename);

% pull up gui to inspect all units in your project
cell_metrics = CellExplorer('metrics',cell_metrics);


%% find ripples
%hpc04 =[38,33] 
%hpc01 = [23,29]
% epochs = [];
% for i = 1:length(session.epochs)
%     epochs = [epochs;session.epochs{i}.startTime session.epochs{i}.stopTime];
% end

% first input [ripple channel, sharp wave channel] using intan channels. 
ripples = DetectSWR([33,3],'basepath',basepath);


% %%                      Local Function Below

function make_xml(basepath,varargin)

p = inputParser;
addParameter(p,'probe_map','A1x64-Poly2-6mm-23s-160.xlsx');
addParameter(p,'lfp_fs',1250,@isnumeric);
parse(p,varargin{:});
probe_map = p.Results.probe_map;
lfp_fs = p.Results.lfp_fs;
[~,basename] = fileparts(basepath);

addpath(genpath(fullfile('external_packages','buzcode')))
if isempty(dir([basepath,filesep,'amplifier.xml'])) % Make xml file
    
    % Check basepath for xml
%     d = dir([basepath,filesep,'**\*.rhd']);
    
    [~, ~, ~, ~,...
    ~, ~, frequency_parameters,~ ]= ...
    read_Intan_RHD2000_file_snlab(basepath);
    
    % pull recording from rhd file
    fs = frequency_parameters.amplifier_sample_rate;
    
    % make or update xml file
    write_xml(basepath,probe_map,fs,lfp_fs)
    
    % rename amplifier.dat to basename.dat
    if ~isempty(dir([basepath,filesep,'amplifier.dat']))
        disp(['renaming amplifer.dat to ',basename,'.dat'])
        % create command
        command = ['rename ',basepath,filesep,'amplifier.dat',' ',basename,'.dat'];
        system(command); % run through system command prompt
    end


else
    disp('XML found. Check channel mapping/bad channels in Neuroscope then go preprocess this session.')
end
rmpath(genpath(fullfile('external_packages','buzcode')))

end

function preprocess_batch(data_folder)
% checks for evidence of preprocessing in data_folder subfolds and runs
% preprocess_session.m for unprocessed sessions


df = compile_sessions(data_folder);

% loop through basepaths in df and process those that don't have evidence of
% processing (in this case chanMap.mat)
for i = 1:length(df.basepath)
    try
    basepath = df.basepath{i}{:};
    basename = basenameFromBasepath(basepath);
    % Check xml/channel mapping: verify channel map, skip bad channels, and save 
%     if ~exist([basepath,filesep,'amplifier.xml'],'file')
%         make_xml(basepath)
%     end
     
        if isempty(dir(fullfile(basepath,[basename,'.lfp'])))
            preprocess_session(basepath,'digitalInputs',false,'check_epochs',false,'kilosort',false,'overwrite',false)
            channel_mapping('basepath',basepath)
        else
            channel_mapping('basepath',basepath)
            continue
        end
    catch
        disp([basepath, 'processing failed. Check all components for preprocessing present'])
    end
end

end

function write_xml(path,map,Fold,FNew)
% Wrapper for buzcode bz_MakeXMLFromProbeMaps with defaults

defaults.NumberOfChannels = 1;
defaults.SampleRate = Fold;
defaults.BitsPerSample = 16;
defaults.VoltageRange = 20;
defaults.Amplification = 1000;
defaults.LfpSampleRate = FNew;
defaults.PointsPerWaveform = 64;
defaults.PeakPointInWaveform = 32;
defaults.FeaturesPerWave = 4;
[~,basename] = fileparts(path);
bz_MakeXMLFromProbeMaps({map},path,basename,1,defaults)
end



% %% ##########################################################################
%
% %   Sorting using Klusta (NOT WORKING ON ANALYSIS COMPUTER LB 09/14/2021)
% %                       
%
% % ##########################################################################
% %% For Klusta
% d   = dir([basepath,filesep,'*.xml']);
% parameters = LoadXml(fullfile(basepath,d(1).name));
%
% % 5. create klusta_folders as well as prb and prm files
% makeProbeMap(basepath)
%
% % 6. Runs klusta on each shank by calling system
% run_klusta(basepath)
%
% %% 7. Spike sort in Phy (Spike sorting kwik files in phy will update the kwik file)
%
