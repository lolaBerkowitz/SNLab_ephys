% get all the folder paths
function tablegen_script (cohortpath)
    basepath=pathGetter(cohortpath);
    output=table();
    for i=1:length(basepath)
        bp=string(basepath(i));
        output=vertcat(output,tablegen(bp));
    end
    csv = 'meta.csv';
    
    % Write the table to a CSV file
    writetable(output, csv);
    
end



function aviFilePaths = pathGetter (cohortpath)
%Gets all the filepaths of folders that contain avi files
%cohortpath: filepath of the main folder
    tempanimals=dir(fullfile(string(cohortpath),'*'));
    animals=fullfile(cohortpath,{tempanimals(~startsWith({tempanimals.name}, '.')).name});
    [~,~,ext]=fileparts(animals);
    aviFilePaths={};
    if ismember('.avi',ext)
        aviFilePaths=[aviFilePaths,cohortpath];
    elseif isfolder(animals)
        for x=animals
            aviFilePaths=[aviFilePaths,pathGetter(x)]; 
        end
    end
end
% 
% %List of the animal folder paths ****** edit here
% cohortpath='/Volumes/sn data server 3/laura_berkowitz/behavior_validation/appps1_cpp/data/cohort3';
% tempanimals = dir(fullfile(cohortpath,'*'));
% animals=fullfile(cohortpath,{tempanimals(~startsWith({tempanimals.name}, '.')).name});
% num=length(animals);
% basepath={};
% %Fix order of files   *****change depending on alphabetical ordering
% order=[3,1,2];
% for i=1:num
%     path = string(animals(i));  % Specify the folder path
%     file_struct = dir(fullfile(path,'*day*'));
%     paths = fullfile(path, {file_struct.name});
%     paths = paths(order);
%     for j=1:length(paths)    
%         basepath{end+1}=paths(j);
%     end
% end
% output=table();
% for i=1:length(basepath)
%     bp=string(basepath(i));
%     output=vertcat(output,tablegen(bp));
% end
% csv = 'meta.csv';
% 
% % Write the table to a CSV file
% writetable(output, csv);