function process_kilosort(basepath,varargin)
    % Function for processing kilosort
    %multishank=1 when multiple shanks and =0 when single shank

    p = inputParser;
    
    addParameter(p,'overwrite',false,@islogical);
    addParameter(p,'multishank',false,@islogical);
    addParameter(p,'ssd_path',"C:\kilo_temp",@ischar);

    parse(p,varargin{:});

    overwrite = p.Results.overwrite;
    multishank = p.Results.multishank;
    ssd_path = p.Results.ssd_path;


    [~,basename] = fileparts(basepath);
    session = loadSession(basepath,basename);
    ssd_folder = fullfile(ssd_path, basename);

    if isempty(dir([ssd_folder+'\Kilosort*'])) || overwrite == true

        % For Kilosort: create channelmap
        if ~isempty(session.channelTags.Bad.channels)
            createChannelMapFile_KSW(basepath,basename,'staggered',session.channelTags.Bad.channels);
        else
             createChannelMapFile_KSW(basepath,basename,'staggered');
        end
        
        % creating a folder on the ssd for chanmap,dat, and xml
        mkdir(ssd_folder);
        
        % Copy chanmap,basename.dat, and xml
        disp('Copying basename.dat, basename.xml, and channelmap to ssd')
        
        disp('Saving dat file to ssd')
        basename=char(basename);
        basepath=char(basepath);
        ssd_folder=char(ssd_folder);
        command = ['robocopy "',basepath,'" ',ssd_folder,' ',basename,'.dat'];
        system(command);
        
        disp('Saving xml to ssd')
        command = ['robocopy "',basepath,'" ',ssd_folder,' ',basename,'.xml'];
        system(command);
        
        disp('Saving channel_map to ssd')
        command = ['robocopy "',basepath,'" ',ssd_folder,' chanMap.mat'];
        system(command);
        
%         run_ks1(basepath,'ssd_folder',ssd_folder)
        % single sort
        kilosortFolder = KiloSortWrapper('SSD_path', ssd_path, ...
            'rejectchannels', session.channelTags.Bad.channels);
    else
        disp("Kilosort folder already exists. Will not overwrite.")
    end

end