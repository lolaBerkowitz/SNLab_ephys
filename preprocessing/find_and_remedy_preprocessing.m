
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

p.parse(varargin{:})
remedy_check = p.Results.remedy_check;

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
get_sessions(status,remedy_check);


end

function [basepaths,f] = get_sessions(status,remedy_check)

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
        basepaths = status.basepath(status.kilosort == 0);
        for i = 1:length(basepaths)
            disp(['Processing : ',basepaths{i}])
            process_tracking(basepaths{i})
        end
        f.a = @create_channelmap;
        f.b = @run_ks1;
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
        basepaths = status.basepath(status.kilosort == 1 & status.phyRez == 1 & status.cell_metrics == 0);
        
        for i = 1:length(basepaths)
            disp(['Processing : ',basepaths{i}])
            basename = basenameFromBasepath(basepaths{i});
            session = loadSession(basepaths{i},basename);
            ProcessCellMetrics('session',session)
        end
        
    case 'ripples'
        basepaths = status.basepath(status.lfp == 1 & status.ripples == 0);
        f = @findRipples;
        
        for i = 1:length(basepaths)
            disp(['Processing : ',basepaths{i}])
            process_tracking(basepaths{i})
        end
end

end
