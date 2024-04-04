% This script postprocesses behavior data only in CellExplorer format.
function process_behavior_batch(data_path,metadata_path,varargin)


% Assumes general behavior file and *maze_coords.csv is in basepath.
p = inputParser;
p.addParameter('overwrite',false,@islogical)

p.parse(varargin{:});
overwrite = p.Results.overwrite;


% loop through basepaths and run main function
df = compile_sessions(data_path);

for i = 1:length(df.basepath)
    df_sub = compile_sessions(df.basepath{i}{1});
    
    for ii = 1:length(df_sub.basepath)
        basepath = df_sub.basepath{ii}{1};
        cd(basepath)
        basename = basenameFromBasepath(basepath);
        disp(['processing basepath: ',basepath])
        
        if check_folders(basepath)
            if exist(fullfile(basepath,[basename,'.animal.behavior.mat'])) & ~overwrite
                continue
            else
                main(basepath,metadata_path)
            end
        else
            disp('Missing DLC for a video')
            continue
        end
    end
    
end

end


function main(basepath,metadata_path)

basename = basenameFromBasepath(basepath);
% check is session file exists, if not make one 


if  ~exist([basepath,filesep,[basename,'.session.mat']],'file')
    session = sessionTemplate(basepath); 
    save(fullfile(basepath,[basename,'.session.mat']),'session')
end

% Make events from DLC
make_events('basepath',basepath);

% update epochs from digitalIn.events.mat
update_epochs('basepath',basepath,...
    'annotate',true,...
    'overwrite',true,...
    'ttl_method',[])

% general behavior file
general_behavior_file_SNlab('basepath',basepath,'force_overwrite',true,'smooth_factor',.2,'primary_coords_dlc',2:3);

% update behavior file from metadata csv
update_behavior_from_metadata(metadata_path,'basepath',basepath);

get_maze_XY('basepath',basepath,'config_path', 'C:\Users\schafferlab\github\SNLab_ephys\behavior\behavior_configs\');

restrict_and_transform(basepath,'overwrite',true);

end

function check = check_folders(basepath)
% see if DLC files are present for each video

% find all .avi videos and dlc files in basepath
video_files = dir(fullfile(basepath,'*.avi'));
video_files = video_files(~cellfun(@(x) ismember(x(1,2),'._'), {video_files.name}));
video_files = extractBefore({video_files.name},'.avi');

dlc_files = get_dlc_files_in_basepath(basepath);
dlc_files = extractBefore({dlc_files.name},'DLC');

% see if each video has a DLC file
if sum(contains(video_files,dlc_files)) == length(video_files)
    check = true;
else
    check = false;
end


end