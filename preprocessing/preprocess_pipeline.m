% Procedures for processing electrophysiology data in the SNLab
%
% Dependencies
%   - external packages: neurocode,KS2Wrapper
%   - schafferlab/github/Kilosort2.5
%   - schafferlab/github/CellExplorer
%
% Data organization assumptions
%   - Data should be organized with each recordding session saved within an animal folder i.e. 'animal_id/basename'
%   - 'basename' is used by CellExplorer and other buzcode functions.
%     Therefore, it is important to make this name unique. 
%   
% LBerkowitz 2021

% Set paths for data
data_folder = 'D:\spike_sorting_practice\2021-09-16_test_210916_135552';

% Set basename from folder
[~,basename] = fileparts(data_folder);
% probe_map = 'CED_E1_4X16_front_front.xlsx';
ssd_path = 'C:\Kilo_temp';

%% ##########################################################################

%                       Setup data for clustering

% ##########################################################################

%% Preprocessing (make metadata file, filter/CAR before kilosort, run kilosort, analyze in phy, ) .

% 1. creates csv for anatomical position and updates basename.session.info
channel_mapping('basepath',data_folder,'show_gui_session',true,'fig',true)

% 2. get header file info and load
d = dir([data_folder,filesep,'**\*.rhd']);
intanRec = read_Intan_RHD2000_file_snlab([d(1).folder,filesep],d(1).name);

%% 3. Make or update xml file
write_xml(data_folder,'A1x64-Poly2-6mm-23s-160.xlsx',30000,1250)

channel_mapping('basepath',data_folder,'show_gui_session',true,'fig',true)

%% 4. In Neuroscope, check channel map, skip bad channels, and save.
% follow steps for chosen spike sorting method (Kilosort2.5 or Klusta)

%% ##########################################################################

%                       Sorting using Kilosort2.5

% ##########################################################################

%% For Kilosort: Create lfp file and make channel map

% 5. Update channel map from basename.xml
create_channelmap(data_folder)

%% Spike sorting
% Move chanMap, xml, and dat to SSD folder.
%creating a folder on the ssd for chanmap,dat, and xml
ks2_folder = fullfile(ssd_path, basename);
mkdir(ks2_folder);

% 6. Spike sort using kilosort 2.5 (data on ssd)
run_ks2(data_folder,ks2_folder)

% 7. Clean up kilo results in Phy

%% 8. Copy back over to data file
disp('Saving data back to data folder from ssd')
command = ['robocopy ',ks2_folder,' ',data_folder,' /e'];
system(command);


%% 9- extract spike times and waveforms for sorted clusters
session = sessionTemplate(data_folder,'showGUI',false);
f = dir([data_folder,filesep,'Kilosort*']);
spikes = loadSpikes('session',session,'clusteringpath',[f.folder filesep f.name]);

%% 10 - compute basic cell metrics
cell_metrics = ProcessCellMetrics('session',session,'manualAdjustMonoSyn',false);

% GUI to manually curate cell classification. This is not neccesary at this point.
% It is more useful when you have multiple sessions
cell_metrics = CellExplorer('metrics',cell_metrics);

%% 11 - copy data to main session folder
% After you have run all this, you need to copy the Kilosort folder
% (actully only some files are neccesary), plus spikes.cellinfo.mat and cell_metrics.cellinfo.mat
% to the main session folder in the shared network drive

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
% %                       Sorting using Klusta (NOT WORKING ON ANALYSIS COMPUTER
% %                       YET LB 09/14/2021)
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

