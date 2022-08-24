function batch_kilosort(data_folder)

% Find all subfolders 
% handle input as csv of basepaths or directory
if contains(data_folder,'.csv')
    df = readtable(data_folder);
else
    df = compile_sessions(data_folder);
end

basepaths = unique(df.basepath);
% loop through folders and process those that don't have evidence of
% processing (in this case chanMap.mat)
for i = 1:length(basepaths)
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    
    cd(basepath)
    ks_check = dir(fullfile(basepath,'Kilosort_*'));
    xml_check = fullfile(basepath,[basename,'.xml']);
    
    if ~isempty(ks_check)
        disp('kilosort already run, moving to next folder')
        continue
    end
    
    if ~exist(xml_check,'file')
        disp('No xml found. Create and check bad channels')
        continue
    end
    
    % create channelmap 
    create_channelmap(basepath)
    
    % run kilosort save ks folder to basepath
    run_ks1(basepath,'ssd_folder',basepath)
    
end
end