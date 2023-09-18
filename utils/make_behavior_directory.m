function make_behavior_directory(basepaths,task_name)
% make_dataset_directory_from_files takes basepaths containing basepath
% variable with paths to each subject folder and creates subdirectories for
% a specified behavior task. 

% LB 06/30/2023

% Input: 
%    basepaths: path to cohort metadata csv containing paths to each subject folder
% 
% Saves directories for each behavior folder 

% read basepaths 
df = readtable(basepaths,"Delimiter",'comma');

for i = 1:length(df.basepath)
    
    % iterate within basepath 
    subject_path = df.basepath{i};
    
    switch task_name
        case 'cpp'
            task_phases = {'pretest_pairing_day01','pairing_day02','pairing_posttest_day03'};
        case 'object-location'
            task_phases = {'hab_day01','hab_day02','hab_day03','object_location_day04'};
        case 'novel-object'
            task_phases = {'hab_day01','hab_day02','hab_day03','novel_object_day04'};
        case 'y-maze'
            task_phases = {'y_maze_alternation'};
        case 't-maze'
            task_phases = {'hab_day01','delay_day01','delay_day02','delay_day03','delay_day04','delay_day05','delay_day06'};
        case 'w-maze'
            task_phases = {'hab_day01','single_arm_day01','single_arm_day02','single_arm_day03',...
                'alternation_day01','alternation_day02','alternation_day03','alternation_day04',...
                'alternation_day05','alternation_day06','alternation_day07'};
    end
    
    % Makes basepaths from names in task_phases
    make_directories(subject_path,task_phases)
    
    % Makes a basename.session file for each basepath in subject_path
    initialize_dir_for_CellExplorer(subject_path)
end


end

function make_directories(subject_path,task_phases)

% create directories for each phase of the task and place in basepath
basename = basenameFromBasepath(subject_path);

for i = 1:length(task_phases)
    % make a directory named basename_task_phase for each item in
    % task_phases variable
    if ~isfolder(fullfile(subject_path,[basename,'_',task_phases{i}]))
        mkdir(fullfile(subject_path,[basename,'_',task_phases{i}]))
    end
end

end

function initialize_dir_for_CellExplorer(subject_path)

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