% get all the folder paths
num=6; %num of animals
%List of the animal folder paths ****** edit here
animals={'/Volumes/sn data server 3/laura_berkowitz/behavior_validation/appps1_cpp/data/cohort2/3768',...
    '/Volumes/sn data server 3/laura_berkowitz/behavior_validation/appps1_cpp/data/cohort2/3769',...
    '/Volumes/sn data server 3/laura_berkowitz/behavior_validation/appps1_cpp/data/cohort2/3770',...
    '/Volumes/sn data server 3/laura_berkowitz/behavior_validation/appps1_cpp/data/cohort2/3775',...
    '/Volumes/sn data server 3/laura_berkowitz/behavior_validation/appps1_cpp/data/cohort2/3778',...
    '/Volumes/sn data server 3/laura_berkowitz/behavior_validation/appps1_cpp/data/cohort2/3779'...
    };
basepath={};
%Fix order of files   *****change depending on alphabetical ordering
order=[3,1,2];
for i=1:num
    path = string(animals(i));  % Specify the folder path
    file_struct = dir(fullfile(path,'*day*'));
    paths = fullfile(path, {file_struct.name});
    paths = paths(order);
    for j=1:length(paths)    
        basepath{end+1}=paths(j);
    end
end
output=table();
for i=1:length(basepath)
    bp=string(basepath(i));
    output=vertcat(output,tablegen(bp));
end
csv = 'meta.csv';

% Write the table to a CSV file
writetable(output, csv);