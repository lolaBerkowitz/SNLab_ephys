function initialize_dir_for_CellExplorer(subject_path)
%% For a subject directory, this function will initalize all basepaths within the subject folder with basename.session.mat files. 

% Use dir to get sessions in subject folder 
 session_list = compile_sessions(subject_path);

% loop through session and create a basename.session file
for i = 1:length(session_list.basepath)
    
    basepath = session_list.basepath{i}; 
    basename = session_list.basename{i};
    
    if ~exist(fullfile(basepath,[basename, '.session.mat']),'file')
        session = sessionTemplate(basepath,'basename',basename);
        save(fullfile(basepath,[basename, '.session.mat']),'session');
    else
        disp('session file already created, skipping basepath')
    end
    
end


end