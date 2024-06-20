function process_kilosort(basepath,multishank)
    % Function for processing kilosort
    %multishank=1 when multiple shanks and =0 when single shank
    ssd_path = 'C:\kilo_temp';
    [~,basename] = fileparts(basepath);
    session = loadSession(basepath,basename);

    if ~multishank
        % for single poly2 shank 64channels
        create_channelmap(basepath)
    
    elseif multishank
        % For Kilosort: create channelmap
        if ~isempty(session.channelTags.Bad.channels)
            createChannelMapFile_KSW(basepath,basename,'staggered',session.channelTags.Bad.channels);
        else
             createChannelMapFile_KSW(basepath,basename,'staggered');
        end
    end
    
    % creating a folder on the ssd for chanmap,dat, and xml
    ssd_folder = fullfile(ssd_path, basename);
    mkdir(ssd_folder);
    
    % Copy chanmap,basename.dat, and xml
    disp('Copying basename.dat, basename.xml, and channelmap to ssd')
    
    disp('Saving dat file to ssd')
    command = ['robocopy "',basepath,'" ',ssd_folder,' ',basename,'.dat'];
    system(command);
    
    disp('Saving xml to ssd')
    command = ['robocopy "',basepath,'" ',ssd_folder,' ',basename,'.xml'];
    system(command);
    
    disp('Saving channel_map to ssd')
    command = ['robocopy "',basepath,'" ',ssd_folder,' chanMap.mat'];
    system(command);
    
    % Spike sort using kilosort 1 (data on ssd)
    run_ks1(basepath,'ssd_folder',ssd_folder)
end