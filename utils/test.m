data_folder = 'E:\data\hpc04';
folders = dir(data_folder);
folders = folders(~ismember({folders.name},{'.','..'}),:);

% loop through folders and process those that don't have evidence of
% processing (in this case chanMap.mat)
for i = 1:length(folders)
    basepath = [data_folder,filesep,folders(i).name];
    basename = basenameFromBasepath(basepath);
    cd(basepath)
%     check = dir([basename,'.ripples.events.mat']);
%     
%     if ~isempty(check)
%         disp('ripples already created')
%         continue
%     end
        
    if ~isempty(dir(fullfile(basepath,[folders(i).name,'.lfp']))) && ...
            ~isempty(dir(fullfile(basepath,'anatomical_map.csv')))
        
        % updates basename.session with channel map
        channel_mapping('basepath',basepath,'fig',false)
        session = loadSession(basepath,basename);
        
        try
            ripple_channel = session.brainRegions.CA1sp.channels(end-1);
            
            ripples = FindRipples('basepath',basepath,...
                'channel',ripple_channel,...
                'thresholds',[0.5 1.5],...
                'durations',[50 500]);
            
            save(fullfile(basepath,[basename '.ripples.events.mat']),'ripples')
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
        end
        
    end
end