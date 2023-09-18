subID = 'hpc07';
data_folder = fullfile('Y:\laura_berkowitz\app_ps1_ephys\data\',subID);
metadata_path = 'Y:\laura_berkowitz\app_ps1_ephys\behavior\behavior_metadata.csv';

df = compile_sessions(data_folder);
meta = readtable(metadata_path,'Delimiter','comma');

meta = meta(contains(meta.basepath,subID),:);
[basepaths,a,c] = unique([meta.basepath]);
% loop through folders and process those that don't have evidence of
% processing (in this case chanMap.mat)

%% %% Check epochs for behavior sessions and update paradigm
context_paths = readtable("Y:\laura_berkowitz\app_ps1_ephys\analysis\context_sessions.csv",'Delimiter','comma');
basepaths = context_paths.basepath;

meta_subset = meta(ismember(meta.basename,basenameFromBasepath(basepaths)),:);
for i = 20:length(basepaths)
%     basepath = df.basepath{i};
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    load(fullfile(basepath,[basename,'.animal.behavior.mat']))
    
    idx = find(ismember(meta_subset.basename,basename));

    temp_paradigm = meta_subset.paradigm{idx};
    for id = 1:length(temp_paradigm)
        behavior.trialsID{id,1} = temp_paradigm(id);
    end
%    gui_session(basepath)
   
   
   save(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')

end

    



%% Check epochs for behavior sessions and update paradigm

for i = 1:length(basepaths)
        
%     basepath = df.basepath{i};
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
% load animal behavior file 
try
    load(fullfile(basepath,[basename,'.animal.behavior.mat']))
catch
    disp('No behavior file found. Generating one now if DLC exists')
   if ~isempty(dir(fullfile(basepath,'*DLC*.csv'))) 
           gui_session(basepath)

       % processing the tracking and generates behavior file
       process_tracking(basepath)
      
       % update behavior file with metadata
       update_behavior_from_metadata(metadata_path,'basepath',basepath);
   end
end

if exist('behavior','var')
    disp('updating metadata')
    update_behavior_from_metadata(metadata_path,'basepath',basepath);
end
    

end

%% Check for tracking

for i = 1:length(basepaths)
%     basepath = df.basepath{i};
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
% load animal behavior file 
try
    load(fullfile(basepath,[basename,'.animal.behavior.mat']))
catch
    disp('No behavior file found. Generating one now if DLC exists')
   if ~isempty(dir(fullfile(basepath,'*DLC*.csv'))) 
           gui_session(basepath)

       % processing the tracking and generates behavior file
       process_tracking(basepath)
      
       % update behavior file with metadata
       update_behavior_from_metadata(metadata_path,'basepath',basepath);
   end
end

if exist('behavior','var')
    disp('updating metadata')
    update_behavior_from_metadata(metadata_path,'basepath',basepath);
end
    

end
