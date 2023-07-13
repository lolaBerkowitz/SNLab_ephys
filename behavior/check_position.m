metadata_path = 'Y:\laura_berkowitz\app_ps1_ephys\behavior\behavior_metadata.csv';

meta = readtable(metadata_path,"Delimiter",',');
basepaths = unique(meta.basepath);
basepaths = basepaths(~contains(basepaths,'hpc01')); % remove for now
basenames = basenameFromBasepath(basepaths);

for i = 1:length(basepaths)
    
    % check date of animal behavior file 
    % if Ran in 2022, overwrite
    basepath = basepaths{i};
    basename = basenames{i};
    session = loadSession(basepath,basename);
    behav = dir(fullfile(basepath,[basename,'.animal.behavior.mat']));
    
    % skip basepaths that don't have a behavior file or were created after
    % the edits made in 2023
    if isempty(behav)
        process_tracking(basepath)
        continue
    end
    
    % run check - will overwrite if position coordinates are not scaled to cm 
    % or centered at 0,0
    check_position_scaling(basepath)
    load(fullfile(basepath,[basename,'.animal.behavior.mat']))
    tracking.restrict(session,behavior,basepath);
   
end