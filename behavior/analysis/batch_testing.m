data_folder = 'Y:\laura_berkowitz\alz_stim\data\crimp';
folders = dir(data_folder);
folders = folders(~ismember({folders.name},{'.','..'}),:);
folders = folders(~ismember({folders.name},{'session_check.csv','._session_check.csv'}),:);


% loop through folders and process those that don't have evidence of
% processing (in this case chanMap.mat)d
for i = 1:length(folders)
    basepath = [data_folder,filesep,folders(i).name];
    basename = basenameFromBasepath(basepath);
    cd(basepath)
    
    if ~exist(fullfile(basepath,[basename,'.animal.behavior.mat']),'file')
        continue    
    end
    
    load(fullfile(basepath,[basename,'.animal.behavior.mat']))
    
    if  range(behavior.speed) > 150 
        [velocity, acceleration,distance_vector] = linear_motion(behavior.position.x,behavior.position.y,behavior.sr,.1);
        behavior.speed = velocity; 
        behavior.accelration = acceleration; 
        
        save(fullfile(basepath,[basename,'.animal.behavior.mat']))
    else 
        continue
    end
        
end