
function find_and_remedy_preprocessing(data_folder,varargin)
% find_and_remedy_preprocessing uses session_status to preprocess sessions with missing
% analyses. Analyses carried out relative to completion of dependencies.
% For instance, sleep states can only be run if the .lfp file is present in
% the data bath.
%
% input:
%   data_folder: path to parent directory containing basepaths or csv of
%   basepaths with column basepath.

p = inputParser;
p.addParameter('remedy_check','all')
addParameter(p,'ssd_path','C:\kilo_temp',@ischar);    % Path to SSD disk.

p.parse(varargin{:})
remedy_check = p.Results.remedy_check;
ssd_path = p.Results.ssd_path;

% handle input as csv of basepaths or directory
if contains(data_folder,'.csv')
    df = readtable(data_folder);
else
    df = compile_sessions(data_folder);
end


% first run a check
check_processing_status(data_folder)
data_folder = fileparts(data_folder);

% load check
check = dir(fullfile(data_folder,'session_check*.csv'));
status = readtable(fullfile(data_folder,check.name));

%% remedy processing
get_sessions(status,remedy_check,ssd_path);


end

function [basepaths,f] = get_sessions(status,remedy_check,ssd_path)

switch remedy_check
    case 'tracking'
        basepaths = status.basepath(status.tracking_dlc == 1 & (status.tracking_animalBehavior == 0 | ...
            status.tracking_restrictxy == 0 | ...
            status.tracking_sessionBehavioralTracking == 0));
        for i = 1:length(basepaths)
            disp(['Processing : ',basepaths{i}])
            process_tracking(basepaths{i})
        end
    case 'kilosort'
        basepaths = status.basepath(status.sorting_Kilosort == 0);
        for i = 1:length(basepaths)
            disp(['Processing : ',basepaths{i}])
            basepath = basepaths{i};
            basename = basenameFromBasepath(basepath);
            % create channelmap
            create_channelmap(basepath)
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
            
            % run kilosort save ks folder to basepath
            run_ks1(basepath,'ssd_folder',ssd_folder)
        end
        
    case 'sleep_states'
        basepaths = status.basepath(status.sleep_states == 0 & status.lfp == 1);
        
        for i = 1:length(basepaths)
            disp(['Processing : ',basepaths{i}])
            process_tracking(basepaths{i})
        end
        
        f.a = @SleepScoreMaster; % takes lfp in base 0
        f.b = @thetaEpochs;
    case 'lfp'
        basepaths = status.basepath(status.lfp == 0);
        f = @preprocess_session; % if no lfp, preprocessing has likely not been run.
        
        for i = 1:length(basepaths)
            disp(['Processing : ',basepaths{i}])
            process_tracking(basepaths{i})
        end
    case 'cell_metrics'
        basepaths = status.basepath(status.sorting_Kilosort == 1 & status.sorting_phyRez == 1 & status.cell_metrics == 0);
        
        for i = 1:length(basepaths)
            disp(['Processing : ',basepaths{i}])
            session = sessionTemplate(basepaths{i});
            ProcessCellMetrics('session',session,'basepath',basepaths{i})
            close all
        end
        
    case 'ripples'
        basepaths = status.basepath(status.lfp == 1 & status.ripples == 0);
        
        for i = 1:length(basepaths)
            disp(['Processing : ',basepaths{i}])
            process_tracking(basepaths{i})
        end
end

end
