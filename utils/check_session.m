data_folder = 'Y:\laura_berkowitz\app_ps1_ephys\data\hpc06';
folders = dir(data_folder);
folders = folders(~ismember({folders.name},{'.','..'}),:);

% loop through folders and process those that don't have evidence of
% processing (in this case chanMap.mat)
for i = 1:length(folders)
    basepath = [data_folder,filesep,folders(i).name];
    basename = basenameFromBasepath(basepath);
    cd(basepath)
    
    session = loadSession(basepath,basename);
%     gui_session
    cell_metrics = CellExplorer('basepath',basepath) 
    
end