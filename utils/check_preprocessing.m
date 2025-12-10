df = compile_sessions("Y:\laura_berkowitz\behavior_metadata.csv");

rewrite_bool = zeros(length(df.basepath),1);

for i = 1:length(df.basepath)

    basepath = df.basepath{i};

    % load chanMap
    if exist(fullfile(basepath,'chanMap.mat'),'file')
        
        load(fullfile(basepath,'chanMap.mat'))
        
        if any(ismember(kcoords,[2,3,4,5]))
            rewrite_bool(i) = true;
        end
    end
end


