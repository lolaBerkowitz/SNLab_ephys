function preprocess_batch_from_df(sessions_csv_path)

sessions = compile_sessions(sessions_csv_path);

for i = 1:length(sessions.basepath)

    basepath = sessions.basepath{i};

    if ismember(sessions.shank_type{i},'single')

        preprocess_session(basepath,'digitalInputs',true,'kilosort',false,'tracking',false,'check_epochs',false)

    else
        preprocess_session(basepath,'digitalInputs',false,'kilosort',false,'tracking',false,'specialChannels',[],'multishank',true,'check_epochs',false)
    end

end