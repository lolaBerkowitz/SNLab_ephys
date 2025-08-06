function make_behavior_directory(subid_csv,data_dir,task_name)
% make_dataset_directory_from_files takes basepaths containing basepath
% variable with paths to each subject folder and creates subdirectories for
% a specified behavior task. 

% LB 06/30/2023

% Input: 
%    subid_csv (String): csv with subject id as a column variable named 'subid' 
%    data_dir (Char): path to main dataset where subject folders can be
%    stored
%    task_name (Char): name of behavioral experiment 
%    basepaths: path to cohort metadata csv containing paths to each subject folder
% 
% Saves directories for each behavior folder 
% subid_csv = 'Y:\laura_berkowitz\behavior_validation\appps1_cheeseboard\subject_list.csv'
% 
% data_dir = 'Y:\laura_berkowitz\behavior_validation\appps1_cheeseboard\data'
% 
% task_name = 'cheeseboard'
% 
% make_behavior_directory(subid_csv,data_dir,task_name)
%

% read basepaths 
df = readtable(subid_csv,"Delimiter",'comma');


for i = 1:length(df.subid)
    
    % iterate within basepath 
    if isa(df.subid(i),'cell')
        subid = df.subid{i};
    elseif isa(df.subid(i),'double')
        subid = num2str(df.subid(i));
    end
    
    subject_path = fullfile(data_dir,subid);
    
    switch task_name
        case 'cpp'
            task_phases = {'habituation_day01','habituation_day02','habituation_day03','cpptask_day04'};
        case 'open_field'
            task_phases = {'open_field_day01','open_field_day02','open_field_day03'};
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
        case 'cheeseboard'
            task_phases = {'hab_day01','hab_day02','hab_day03','reward_day03','reward_day04','reward_day05','task_day06'};
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

