data_folder = 'Y:\laura_berkowitz\app_ps1_ephys\data\hpc05';
folders = dir(data_folder);
folders = folders(~ismember({folders.name},{'.','..'}),:);

% loop through folders and process those that don't have evidence of
% processing (in this case chanMap.mat)
for i = 1:length(folders)
    basepath = [data_folder,filesep,folders(i).name];
    basename = basenameFromBasepath(basepath);
    cd(basepath)
    if ~isempty(dir(fullfile(basepath,[folders(i).name,'.lfp']))) && ...
            ~isempty(dir(fullfile(basepath,'anatomical_map.csv')))
        
        % updates basename.session with channel map
        channel_mapping('basepath',basepath)
        session = loadSession(basepath,basename);
        
        try
            ripple_channel = session.brainRegions.CA1sp.channels(end);
            ripples = FindRipples('basepath',basepath,'channel',ripple_channel,'thresholds',[0.25 1]);
            save([basename '.ripples.events.mat'],'ripples')
        catch
            disp('CA1sp not found, skipping session')
            continue
        end
        
        % refine ripple detection using spiking level
        if ~isempty(dir(fullfile(basepath,[folders(i).name,'.spikes.cellinfo.mat'])))
            spikes = importSpikes('cellType', "Pyramidal Cell", 'brainRegion', "CA1");
            ripplesTemp = eventSpikingTreshold(ripples,'spikes',spikes,'spikingThreshold',0.5);
            ripples = ripplesTemp;
            if ~exist([basepath '\Ripple_Profile'])
                mkdir('Ripple_Profile');
            end
            save([basename '.ripples.events.mat'],'ripples');saveas(gcf,['Ripple_Profile\SWRmua.png']);
        end
        
    end
end