% load sessions csv 
metadata = readtable('Y:\laura_berkowitz\behavior_metadata.csv','Delimiter','comma');

% only use the basepaths 
basepaths = unique(metadata.basepath);

%% overwrite 
overwrite_coords = true; 

for i = 1:length(basepaths)
    % get channel map from anatomical 
    basepath = basepaths{i};
    cd(basepath)
    channel_mapping('basepath',basepath,'fig',false)
    
end
    [~, basename] = fileparts(basepath);
    
    temp_metadata = metadata(contains(metadata.basepath,basepath),:);

    if ~exist(fullfile(basepath,[basename,'.session.mat']))
        session = sessionTemplate(basepath);
    else
        session = loadSession(basepath,basename);
    end

    temp_metadata = metadata(contains(metadata.basepath,basepath),:);
    probe_type = temp_metadata.probe_type(1);

    if contains(probe_type,'A1x64')

        coords = readtable('C:\Users\schafferlab\github\SNLab_ephys\preprocessing\probe_maps\chanCoordinates_A1x64-Poly2-6mm-23s-160.csv');
        shank_spacing = 0;
        vert_space = 46;
        layout = 'poly 2';

    elseif contains(probe_type,'A5x12')

        coords = readtable('C:\Users\schafferlab\github\SNLab_ephys\preprocessing\probe_maps\chanCoordinates_A5x12-16-Buz-lin-5mm-100-200-160-177.csv');
        shank_spacing = 300;
        vert_space = 26;
        layout = 'poly 2';
    end

    % save coords to session file 
    session.extracellular.chanCoords.x = coords.x_um;
    session.extracellular.chanCoords.y = coords.y_um;
    session.extracellular.chanCoords.layout = layout;
    session.extracellular.chanCoords.verticalSpacing = vert_space;
    save(fullfile(basepath,[basename,'.session.mat']),'session')

    % check for 
    try     
        load(fullfile(basepath,[basename,'.chanCoords.channelInfo.mat']))
    catch 
        disp('No chanCoords found. Writing file now')
    end

    if overwrite_coords
    chanCoords.x = coords.x_um;
    chanCoords.y = coords.y_um;
    chanCoords.source = 'SNLab_probe_maps';
    chanCoords.layout = layout;
    chanCoords.shankSpacing = shank_spacing;
    chanCoords.verticalSpacing = vert_space;
    

    % save both back to basepath 
    save(fullfile(basepath,[basename,'.chanCoords.channelInfo.mat']),'chanCoords')
    end
end
