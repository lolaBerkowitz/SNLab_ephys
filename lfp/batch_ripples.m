function batch_ripples(data_folder,varargin)
% For all session files in Data folder (path to animal folder or csv of basepaths), runs
% findRipple algorithm with default paramters. If session has cell_metrics,
% will also run eventSpikingThreshold to keep keep ripples that have
% increased firing of CA1 units greater than .5. 

p = inputParser;
p.addParameter('overwrite',false,@islgical)
p.addParameter('thresholds',[0.1 .5])
p.addParameter('durations',[50 500])
p.addParameter('minDuration',20)
p.addParameter('spikingThreshold',0.5)

p.parse(varargin{:});
overwrite = p.Results.overwrite;
thresholds = p.Results.thresholds;
durations = p.Results.durations;
minDuration = p.Results.minDuration;
spikingThreshold = p.Results.spikingThreshold;

% if the data_folder is a path leading to csv load with readtable
if contains(data_folder,'.csv')
    df = readtable(data_folder);
    [sessions.folder,sessions.name] = fileparts(df.basepath);
else % load with dir
    sessions = dir(data_folder);
    sessions = sessions(~ismember({sessions.name},{'.','..'}),:);
    sessions = sessions([sessions.isdir]);
end

% loop through folders and process those that don't have evidence of
% processing (in this case chanMap.mat)
for i = 1:length({sessions.name})
%     basepath = df.basepath{i};
    basepath = fullfile(sessions(i).folder,sessions(i).name);
    basename = sessions(i).name;
    
    % updates basename.session with channel map
    channel_mapping('basepath',basepath,'fig',false)
       
    cd(basepath) % do this because certain functions have insidious pwd scattered throught code
    
    % skip session if overwrite is false and ripples have already been
    % created
    if exist(fullfile(basepath,[basename '.ripples.events.mat']),'file')
        if overwrite
            continue
        end
    end
        
    % only run sessions with lfp file created and 
    if ~isempty(dir(fullfile(basepath,[basename,'.lfp']))) && ...
            ~isempty(dir(fullfile(basepath,'anatomical_map.csv')))
        
        % load session to find CA1 channels
        session = loadSession(basepath,basename);

        try
            ripple_channel = session.brainRegions.CA1sp.channels(end-1);
            ripples = FindRipples('basepath',basepath,...
                'channel',ripple_channel,...
                'thresholds',thresholds,...
                'durations',durations,...
                'minDuration',minDuration);
        catch
            % above will fail if no CA1 units detected
            disp('CA1sp not found, skipping session')
            continue
        end
        
        % refine ripple detection using CA1 spiking
        if ~isempty(dir(fullfile(basepath,[basename,'.spikes.cellinfo.mat'])))
            spikes = importSpikes('basepath',basepath,'brainRegion', "CA1");
            ripplesTemp = eventSpikingTreshold(ripples,'spikes',spikes,'spikingThreshold',spikingThreshold);
            ripples = ripplesTemp;
            if ~exist([basepath '\Ripple_Profile'])
                mkdir([basepath '\Ripple_Profile']);
            end
            save(fullfile(basepath,[basename '.ripples.events.mat']),'ripples');
            saveas(gcf,fullfile(basepath,['Ripple_Profile\SWRmua.png']));
        end
        
    end
end