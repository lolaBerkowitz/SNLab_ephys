function check_process_summary(paths)

%Run to update summary for all folders titled data
%paths is a list (>=1) of strings of parent folder paths containing target
%"data" folders
%
%Input (for mac):
%   paths (str array): ["/Volumes/sn data server 3/laura_berkowitz/app_ps1_ephys",...
%   "/Volumes/sn data server 3/laura_berkowitz/alz_stim"];
%
%Output:
%   Runs check_processing_status for each mouse in parent folder and
%   creates Combined.csv in folder containing check_process_summary.m which
%   combines all session_check.csv created by check_processing_status

% updated documention to clarify the output csv 

filePaths=[];
for i=1:length(paths)
    %Calls check_processing_status for each folderpath under data folders
    pathList=getPaths(paths(i));
    filePaths=[filePaths,pathList];

    for i=1:length(pathList)
        check_processing_status(char(pathList(i)));
    end
end
%Create the summary table from all session_check csvs in the specified
%paths
summaryTable(filePaths);
%TODO insert path of .bat file and uncomment line below if running python
%code from here
%system(INSERT,'-echo');
end


%Recursive helper function to get all folder pathnames under data folders
function filePaths=getPaths(path)
filePaths={};
fileName=dir(path);
for i=1:length(fileName)
    temp=fileName(i);
    if strcmp(temp.name,'data')
        dataPath=dir(fullfile(path,temp.name));
        for j=1:length(dataPath)
            %Exclude 'to_split' folders
            if dataPath(j).isdir && ~strcmp(dataPath(j).name,'.') && ...
                    ~strcmp(dataPath(j).name,'..') && ~strcmp(dataPath(j).name,'to_split')
                filePaths=[filePaths,fullfile(dataPath(j).folder,dataPath(j).name)];
           end
        end
        return
    elseif temp.isdir && ~strcmp(temp.name,'.') && ~strcmp(temp.name,'..')
        % Recursively call the function for subfolders
        subfolder=fullfile(path,temp.name);
        dataPath=getPaths(subfolder);
        filePaths=[filePaths,dataPath];
    end
end
end


function summaryTable(filePaths)
%Helper function that takes all session_check CSVs and makes summary table excluding dates

dataCell={};
% Read and combine data from all 'session_check.csv' files
for i=1:length(filePaths)
    filePath=fullfile(char(filePaths(i)), 'session_check.csv');
    data=readtable(filePath);
    columnsDel=[4,6,8,10,12,14,16,18,20,22,24,26];
    data(:,columnsDel)=[];
    % Store data in the cell array
    dataCell{end+1} = data;
end

% Combine data tables into one table
combinedData=vertcat(dataCell{:});
scriptFolder=fileparts(mfilename('fullpath'));
%TODO change / to \ for windows
outputFile=[scriptFolder,'\Combined.csv'];
% Write the combined data to a single CSV file
writetable(combinedData,outputFile);
end
