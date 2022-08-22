% batch_preprocessing uses session_status to preprocess sessions with missing
% analyses. Analyses carried out relative to completion of dependencies.
% For instance, sleep states can only be run if the .lfp file is present in
% the data bath.

data = 'Y:\laura_berkowitz\app_ps1_ephys\data';
remedy_check = 'tracking';

folders = dir(data);
folders = folders(~ismember({folders.name},{'.','..'}),:); % remove ., ..
folders = folders(~ismember({folders.name},'to_split'));

for i = 1:length({folders.name})
    data_folder = fullfile(folders(i).folder,folders(i).name);
    
    % first run a check
    check_processing_status(data_folder)
    
    % load check
    check = dir(fullfile(data_folder,'session_check*.csv'));
    status = readtable(fullfile(data_folder,check.name));
    
    % remedy processing 
    [basepaths,func] = get_sessions(status,remedy_check);
    
    for sess = 1:length(basepaths)
        % run function 
        func('basepath',basepath) 
    end
    
end



function [basepaths,f] = get_sessions(status,remedy_check)

switch remedy_check
    case 'tracking'
        basepaths = status.basepath(status.tracking_dlc == 1 & status.tracking_animalBehavior == 0);
        f = @process_tracking;
    case 'kilosort'
        basepaths = status.basepath(status.kilosort == 0);
        f.a = @create_channelmap;
        f.b = @run_ks1;
    case 'sleep_states'
        basepaths = status.basepath(status.sleep_states == 0 & status.lfp == 1);
        f.a = @SleepScoreMaster; % takes lfp in base 0
        f.b = @thetaEpochs;
    case 'lfp'
        basepaths = status.basepath(status.lfp == 0);
        f = @preprocess_session; % if no lfp, preprocessing has likely not been run. 
    case 'cell_metrics'
        basepaths = status.basepath(status.kilosort == 1 & status.phyRez == 1 & status.cell_metrics == 0);
        f = @ProcessCellMetrics;
    case 'ripples'
        basepaths = status.basepath(status.lfp == 1 & status.ripples == 0);
        f = @findRipples;
end

end

