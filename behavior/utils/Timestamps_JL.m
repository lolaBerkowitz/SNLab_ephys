
function outTable=Timestamps_JL(basefolder)
%Function to get the difference between the start times of the videos in
%the folder and each subfolder contained in the main folder
%basefolder: filepath of the main folder. Any .avi videos in the main
%folder will not be in the csv generated

% Get all the videos from the directory
    paths=pathGetter(basefolder);
    outTable=table();
    maxnum=maxFinder(paths);
    % Get table headers
    File_name={'File name'};
    for i=1:maxnum-1
        new="Start "+string(i+1)+" - Start "+i;
        File_name{1+i}=new;
    end
    
    for i=1:1:length(paths);
        basepath=string(paths(i));
        file_struct = dir(fullfile(basepath,'*.avi'));
        % row names
        bCharPath=char(basepath);
        File_name=bCharPath(1+max(strfind(basepath,'/')):end); 
        % Extract the timestamps from the filenames
        timestamps = regexp({file_struct.name}, '-\d+', 'match');
        timestamps = cellfun(@(x) extractAfter(x,'-'),timestamps,'UniformOutput',false);
        timestamps = cellfun(@(x) str2double(x),timestamps,'UniformOutput',false);
        
        % Sort the file list based on the timestamps
        [~, sortedIndices] = sort( [timestamps{:}]);
        
        % convert to matrix so you can apply matrix functions 
        timestamps = cell2mat(timestamps); 
        %convert to just HHMMSS
        timestamps=mod(timestamps,1000000);
        hours=floor(timestamps/10000);
        minutes=floor((timestamps-hours*10000)/100);
        seconds=mod(timestamps,100);
        %new matrix of time of day in seconds
        newtimestamps=hours*3600+60*minutes+seconds;
        %Sort and find the difference
        newtimestamps=sort(newtimestamps,'ascend');
        difference=diff(newtimestamps/60);
        

        row={File_name};
        for i=1:1:maxnum-1;
            if length(difference)>=i
                 if  floor(difference(i)/60)>0
                    time=floor(difference(i)/60)+" hrs "+ round(mod(difference(i),60)) +" mins";
                 else
                    time=round(mod(difference(i),60)) +" mins";
                 end
                row=[row,time];
            else
                row=[row,NaN];
            end
        end
        row=table(row);
        outTable=vertcat(outTable,row);
        
    end

    % Write the table to a CSV file
    csv = 'timestamps.csv';
    %outTable.Properties.VariableNames([1 2 3]) = {'File name','A','B'}
    writetable(dataTable, csv);
end



% get all the folder paths

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

function outnum=maxFinder(paths)
    outnum=0;
    for x=paths
        templen=numel(dir(fullfile(string(x),'*.avi')));

        if templen>outnum
            outnum=templen;
        end
    end
end


        
        

  
    % output=table();
    % for i=1:length(basepath)
    %     bp=string(basepath(i));
    %     output=vertcat(output,tablegen(bp));
    % end
    % csv = 'timestamps.csv';
    % 
    % % Write the table to a CSV file
    % writetable(output, csv);







% load video
% vid_obj = VideoReader('/Volumes/sn data server 3/laura_berkowitz/behavior_validation/appps1_cpp/data/cohort2/3769/3769_pretest_pairing_day01/3769L_pairing1_A-07122023111233.avi');
% fs = vid_obj.FrameRate;
% duration = vid_obj.Duration;


% function out = pathGetter (bigpath)
%     cohortpath='/Volumes/sn data server 3/laura_berkowitz/behavior_validation/appps1_cpp/data/cohort3';
%     tempanimals = dir(fullfile(cohortpath,'*'));
%     animals=fullfile(cohortpath,{tempanimals(~startsWith({tempanimals.name}, '.')).name});
%     num=length(animals);
%     basepath={};
%     %Fix order of files   *****change depending on alphabetical ordering
%     order=[3,1,2];
%     for i=1:num
%         path = string(animals(i));  % Specify the folder path
%         file_struct = dir(fullfile(path,'*day*'));
%         paths = fullfile(path, {file_struct.name});
%         paths = paths(order);
%         for j=1:length(paths)    
%             basepath{end+1}=paths(j);
%         end
%     end
%     out=basepath;
% end



% Given the format of the timestamps is MMDDYYYYHHMMSS, extract the time
% information so you can compute the inter-video interval or time between
% the start time of each video. 