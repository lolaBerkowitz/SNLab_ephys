
df = compile_sessions('Y:\laura_berkowitz\behavior_validation\appps1_cheeseboard\completed_sessions_10072025.csv');

     
for i = 1:length(df.basepath)

    basepath = df.basepath{i};
    basename = basenameFromBasepath(basepath);
    session = loadSession(basepath,basename);

    for ii = 1:length(session.epochs)
        vidname = extractBefore(session.behavioralTracking{1,ii}.notes,'.');
        rewards_table = session.behavioralTracking{1,ii}.cheeseboard_rewards;
        start_table = session.behavioralTracking{1,ii}.cheeseboard_start;

         % save data to csv 
        reward_file = fullfile(basepath,[vidname,'_rewards.csv']);
        writetable(rewards_table,reward_file);

                 % save data to csv 
        start_file = fullfile(basepath,[vidname,'_start.csv']);
        writetable(start_table,start_file);
    end
end

% basepath = 'Y:\laura_berkowitz\behavior_validation\appps1_cheeseboard\data\cohort_2\4441\4441_exposure3_phase3'
% behavior_funcs.scale_napari_coords(basepath,'csv_tag', '_rewards')
% 
% Scale XY coordinates of rewards, save to session and basepath
% behavior_funcs.scale_napari_coords(basepath,'csv_tag','_start')