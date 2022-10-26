data_folder = 'Y:\laura_berkowitz\app_ps1_ephys\data\hpc10';
metadata_path = 'Y:\laura_berkowitz\app_ps1_ephys\behavior\behavior_metadata.csv';

df = compile_sessions(data_folder);
basepaths = unique([df.basepath{:}])';
% loop through folders and process those that don't have evidence of
% processing (in this case chanMap.mat)

% basepaths(contains(basepaths,'hpc01')) = [];
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
