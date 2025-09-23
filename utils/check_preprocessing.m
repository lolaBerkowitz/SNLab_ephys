df = compile_sessions("Y:\laura_berkowitz\app_ps1_ephys\complete_sessions_07_14_2025.csv");

rewrite_bool = zeros(length(df.basepath),1);

for i = 1:length(df.basepath)

    basepath = df.basepath{i};

    complete_bool = check_lfp(basepath);

    if ~complete_bool
       rewrite_bool(i) = 1;
    end

end


