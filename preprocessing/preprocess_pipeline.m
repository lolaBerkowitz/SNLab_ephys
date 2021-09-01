% Procedures for processing electrophysiology data in the SNLab
%
% Dependencies
%   - external packages: neurocode, Kilosort2, KS2Wrapper
%
% LBerkowitz 2021

% Set paths for data 
data_folder = pwd;
[~,basename] = fileparts(data_folder);
probe_map = 'CED_E1_4X16_front_front.xlsx';
ssd_path = 'C:\Kilo_temp';

addpath(genpath('D:\Users\BClarkLab\ephys_tools\external_packages\analysis-tools'))

%% ##########################################################################

%                       Setup data for clustering  

% ##########################################################################

%% Preprocessing (organize dat based on channel map) 
% Channels were not organized during acquisition. This opens and saves the
% dat file so channels are organized based on channel map. 

%.0  Reorganize dat file if your data was not mapped when streamed to disk
reorganize_dat_from_chanMap(data_folder,probe_map)

% 1. Get filtered raw dat signal
filter_raw_dat_from_dat

% 2. Apply CAR to remove hypersynchronous events 
applyCARtoDat(['amplifier','.dat'], 64, data_folder);

% 3. Make xml file
write_xml(data_folder,'channel_map_64.xlsx',30000,1250) 

%% 4. In Neuroscope, check channel map, skip bad channels, and save. 

% follow steps for chosen spike sorting method (Kilosort2.5 or Klusta)

%% ##########################################################################

%                       Sorting using Kilosort2.5

% ##########################################################################

%% For Kilosort: Create lfp file and make channel map

% 5. Update channel map from basename.xml
create_channelmap(basepath,varargin)

%% Spike sorting
% Move chanMap, xml, and dat to SSD folder. 
%creating a folder on the ssd for chanmap,dat, and xml
ks2_folder = fullfile(ssd_path, basename);
mkdir(ks2_folder);

% 6. Spike sort using kilosort 2.5 (data on ssd)
run_ks2([ssd_path,filesep,basename])

% 7. Copy back over to data file 
disp('Saving data back to data folder from ssd')
command = ['robocopy ',ssd_path,filesep,basename,' ',data_folder,' /e'];
system(command);

%% 8. Clean up kilo results in Phy


%% ##########################################################################

%                       Sorting using Klusta 

% ##########################################################################
%% For Klusta
d   = dir([data_folder,filesep,'*.xml']);
parameters = LoadXml(fullfile(data_folder,d(1).name));

% 5. create klusta_folders as well as prb and prm files 
makeProbeMap(data_folder)

% 6. Runs klusta on each shank by calling system
run_klusta(data_folder)

%% 7. Spike sort in Phy (Spike sorting kwik files in phy will update the kwik file)

%% 9. After spike sort cleanup 
after_spikesort_cleanup.handle_kwik(data_folder)

% 10. Posprocess 
postprocess('DLC',1)


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
    

