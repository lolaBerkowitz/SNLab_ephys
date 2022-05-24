data_folder = 'Y:\laura_berkowitz\app_ps1_ephys\data\hpc05';
folders = dir(data_folder);
folders = folders(~ismember({folders.name},{'.','..'}),:);
overwrite = true;
% loop through folders and process those that don't have evidence of
% processing (in this case chanMap.mat)
for i = 1:length(folders)
    basepath = [data_folder,filesep,folders(i).name];
    basename = basenameFromBasepath(basepath);
    cd(basepath)
    check = dir([basename,'.interictal_spikes.events.mat']);
    
    if ~isempty(check) & ~overwrite
        disp('interictal_spikes already created')
        continue
    end
        
    if ~isempty(dir(fullfile(basepath,[folders(i).name,'.lfp']))) && ...
            ~isempty(dir(fullfile(basepath,'anatomical_map.csv')))
        
        % updates basename.session with channel map
         channel_mapping('basepath',basepath,'fig',false)
        session = loadSession(basepath,basename);
        bad_channels = session.channelTags.Bad.channels;
        try
            cortex_channels = session.brainRegions.Cortex.channels;
            cortex_channels(ismember(cortex_channels,bad_channels)) = [];
            channel = cortex_channels(1);
            interictal_spikes = FindRipples('basepath',basepath,...
                'channel',channel,...
                'thresholds',[.25 6],...
                'durations',[10 50],...
                'minDuration',10);
            save([basename '.interictal_spikes.events.mat'],'interictal_spikes')
        catch
            disp('Cortex not found, skipping session')
            continue
        end
        
%         % refine ripple detection using spiking level
%         if ~isempty(dir(fullfile(basepath,[folders(i).name,'.spikes.cellinfo.mat'])))
%             spikes = importSpikes('cellType', "Pyramidal Cell", 'brainRegion', "CA1");
%             ripplesTemp = eventSpikingTreshold(ripples,'spikes',spikes,'spikingThreshold',0.5);
%             ripples = ripplesTemp;
%             if ~exist([basepath '\Ripple_Profile'])
%                 mkdir('Ripple_Profile');
%             end
%             save([basename '.ripples.events.mat'],'ripples');saveas(gcf,['Ripple_Profile\SWRmua.png']);
%         end
        
    end
end