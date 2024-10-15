function save_mat_not_as_v7_3(project_folder)
% save_mat_not_as_v7_3 saves mat files saved as v7.0 for easier handling in
% python (i.e. using scipy's loadmat)

% grab subject folders from data
folders = dir(project_folder);
folders = folders(~ismember({folders.name},{'.','..'}),:);
folders(contains({folders.name},'to_split'),:) = [];

% loop through subjects
for i = 1:length(folders)
    subfolders = dir(fullfile(folders(i).folder,folders(i).name));
    subfolders = subfolders(~ismember({subfolders.name},{'.','..'}),:);
    
    % loop through sessions
    for ii = 1:length(subfolders)
        %define basepath and basename 
        basepath = fullfile(subfolders(ii).folder,subfolders(ii).name);
        basename = basenameFromBasepath(basepath);
        
        %load ripples and events 
        try
        load(fullfile(basepath, [basename,'.EMGFromLFP.LFP.mat']));
        catch
            disp(['Error loading file for: ',basepath, ' skipping for now'])
            continue
        end
        save(fullfile(basepath, [basename,'.EMGFromLFP.LFP.mat']),'EMGFromLFP');
    end
end

end