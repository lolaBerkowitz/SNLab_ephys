function check_position_scaling(basepath)
% Checks SNLab behavior file position data to see if coordinates have been
% scaled and translated to be centered at 0,0. If not, scales coordinates
% realtive to maze center. 

basename = basenameFromBasepath(basepath);

% load session and behavior file 
try
    load(fullfile(basepath,[basename,'.session.mat']))
    load(fullfile(basepath,[basename,'.animal.behavior.mat']))
catch
    disp('No behavior or session file. skipping')
    return
end

% Check if coordinates have been scaled by verifying units are in cm
scale_flag = contains(behavior.position.units ,'cm');

% Check if center of x,y are at zero
centered_flag = (median(min(behavior.position.x):max(behavior.position.x)) < 1 && ... % check for X near 0
    median(min(behavior.position.x):max(behavior.position.x)) > -1) && ...
    (median(min(behavior.position.y):max(behavior.position.y)) < 1 && ... % check for Y near 0
     median(min(behavior.position.y):max(behavior.position.y)) > -1);

% Scale coordinates if not in cm
if ~scale_flag || ~centered_flag
    % redo general_behavior_file and processing
    general_behavior_file_SNlab('basepath',basepath,'force_overwrite',true);
    update_behavior_from_metadata('Y:\laura_berkowitz\app_ps1_ephys\behavior\behavior_metadata.csv','basepath',basepath)
    get_maze_XY('basepath',basepath,'overwrite',true);
    restrict_and_transform(basepath,'overwrite',true);
end

end