function initialize_dir_for_CellExplorer(subject_path)
%% For a subject directory, this function will initalize all basepaths within the subject folder with basename.session.mat files. 

% Use dir to get sessions in subject folder 
session_list = dir(subject_path);
session_list = session_list(~ismember({session_list.name},{'.','..'}),:); 
session_list = session_list(~cellfun(@(x) contains(x,'DS_Store'),{session_list.name})); 
% loop through session and create a basename.session file
for i = 1:length(session_list)
    
    basepath = fullfile(subject_path,session_list(i).name); 
    basename = basenameFromBasepath(basepath);
    session = sessionTemplate(basepath,'basename',basename);
    save(fullfile(basepath,[basename, '.session.mat']),'session');
    
end


end