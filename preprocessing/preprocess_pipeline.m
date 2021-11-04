% Procedures for processing electrophysiology data in the SNLab
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
%
% LBerkowitz 2021

% Set paths for data
data_folder = 'D:\app_ps1\data\hpc01\day05_211101_140116';
% Set basename from folder
[~,basename] = fileparts(data_folder);
ssd_path = 'C:\kilo_temp';

% settings for xml creation (if needed)
% probe_map = 'A1x64-Poly2-6mm-23s-160.xlsx';
% lfp_fs = 1250;

%% ##########################################################################

%                       Setup data for clustering

% ##########################################################################

%% Preprocessing (update or create xml, check in neuroscope (manually), update channel map, run kilosort, annotate in phy)

% 0. rename amplifier.dat to basename.dat
if ~isempty(dir([data_folder,filesep,'amplifier.dat']))
    disp(['renaming amplifer.dat to ',basename,'.dat'])
    % create command
    command = ['rename ',data_folder,filesep,'amplifier.dat',' ',basename,'.dat'];
    system(command); % run through system command prompt
end

% lets also rename the xml if present.
if ~isempty(dir([data_folder,filesep,'amplifier.xml']))
    disp(['renaming amplifer.xml to ',basename,'.xml'])
    % create command
    command = ['rename ',data_folder,filesep,'amplifier.xml',' ',basename,'.xml'];
    system(command); % run through system command prompt
elseif isempty(dir([data_folder,filesep,'*.xml'])) % Make xml file
    
    % Check basepath for xml
    d = dir([data_folder,filesep,'**\*.rhd']);
    intanRec = read_Intan_RHD2000_file_snlab([d(1).folder,filesep],d(1).name);
    
    % pull recording from rhd file
    fs = intanRec.frequency_parameters.amplifier_sample_rate;
    
    % make or update xml file
    write_xml(data_folder,probe_map,fs,lfp_fs)
end


%% 1. In Neuroscope, verify channel map, skip bad channels, and save.
% follow steps for chosen spike sorting method (Kilosort)


%% ##########################################################################

%                       Sorting using Kilosort

% ##########################################################################

%% For Kilosort:

% 2. Update channel map from basename.xml
create_channelmap(data_folder)

%% Spike sorting

% 3. creating a folder on the ssd for chanmap,dat, and xml
ssd_folder = fullfile(ssd_path, basename);
mkdir(ssd_folder);

%% 4. Copy chanmap,basename.dat, and xml
disp('Copying basename.dat, basename.xml, and channelmap to ssd')

disp('Saving dat file to ssd')
command = ['robocopy ',data_folder,' ',ssd_folder,' ',basename,'.dat'];
system(command);

disp('Saving xml to ssd')
command = ['robocopy ',data_folder,' ',ssd_folder,' ',basename,'.xml'];
system(command);

disp('Saving channel_map to ssd')
command = ['robocopy ',data_folder,' ',ssd_folder,' chanMap.mat'];
system(command);

%% 5. Spike sort using kilosort 1 (data on ssd)
run_ks1(data_folder,ssd_folder)

%% 7. Clean up kilo results in Phy
% In anaconda prompt, cd to kilosort folder.
% cd ks2_folder i.e. cd C:\kilo_temp\2021-09-16_test_210916_135552
% conda activate phy2
% phy template-gui params.py


%% 8. Copy back over to data file
% Copy the kilosort folder from the ssd_path to the main data folder.
%

%% 9- extract spike times and waveforms for sorted clusters
session = sessionTemplate(data_folder,'showGUI',true);

%% 10 - compute basic cell metrics
cell_metrics = ProcessCellMetrics('session',session,'manualAdjustMonoSyn',false);

% GUI to manually curate cell classification. This is not neccesary at this point.
% It is more useful when you have multiple sessions
cell_metrics = CellExplorer('metrics',cell_metrics);


%% 11 - create downsampled lfp

LFPfromDat(data_folder,'outFs',1250,'useGPU',false);

%% To load multiple sessions into CellExplorer
% After you have multiple sessions spike sorted and processed, you can open
% them all in cell explorer.

% assuming your current directory is (../project/animal/session), this will
% cd to project level (for example: A:\Data\GirardeauG)
cd('..\..')

% make this current directory into variable
data_path = pwd;

% look for all the cell_metrics.cellinfo.mat files
files = dir([data_path,'\**\*.cell_metrics.cellinfo.mat']);

% pull out basepaths and basenames
for i = 1:length(files)
    basepath{i} = files(i).folder;
    basename{i} = basenameFromBasepath(files(i).folder);
end

% load all cell metrics
cell_metrics = loadCellMetricsBatch('basepaths',basepath,'basenames',basename);

% pull up gui to inspect all units in your project
cell_metrics = CellExplorer('metrics',cell_metrics);


%% Local functions below

% write xml
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
% d   = dir([data_folder,filesep,'*.xml']);
% parameters = LoadXml(fullfile(data_folder,d(1).name));
%
% % 5. create klusta_folders as well as prb and prm files
% makeProbeMap(data_folder)
%
% % 6. Runs klusta on each shank by calling system
% run_klusta(data_folder)
%
% %% 7. Spike sort in Phy (Spike sorting kwik files in phy will update the kwik file)
%

