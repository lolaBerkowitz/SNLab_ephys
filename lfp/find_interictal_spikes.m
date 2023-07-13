
% Detection heuristics 
%  * Cortical or non-hippocampal channel (given IEDs filtered at
%  physiological ripple frequency (80 - 250Hz) exhibit high power in this
%  band. Thus non-hippocampal channel are used for first detection. 
%  * Passband filtered at 60 - 80 Hz (per Gelinas et al 2016) 
%  * Durations of events from 10 - 60 ms, events within 10ms of each other
%  are merged
%  * Events at least greater than 3 SD above normalized squared signal

% data_folder = 'Y:\laura_berkowitz\app_ps1_ephys\data\hpc05';
% folders = dir(data_folder);
% folders = folders(~ismember({folders.name},{'.','..'}),:);
% overwrite = true;


function find_interictal_spikes(data_path,varargin)

p = inputParser;
addParameter(p,'basepath',[],@isstr)
addParameter(p,'overwrite',true,@islogical)
addParameter(p,'annotate_only',false,@islogical)


parse(p,varargin{:})
overwrite = p.Results.overwrite;
annotate_only = p.Results.annotate_only;
basepath = p.Results.basepath;


% Load sessions 
if ~isempty(basepath)
    sessions = {basepath};
else
    df = compile_sessions(data_path);
    sessions = [df.basepath{:}];
end
% loop through folders and process those that don't have evidence of
% processing (in this case chanMap.mat)
for i = 1:length(sessions)
   
    basepath = sessions{i};
    basename = basenameFromBasepath(basepath);
        
    % see if interictal_spikes events have been processed
    check = dir(fullfile(basepath,[basename,'.interictal_spikes.events.mat']));
    
    if annotate_only
        
        NeuroScope2('basepath',basepath)
        continue
    end
    
    if ~isempty(check) && ~overwrite
        disp('interictal_spikes already created')
        continue
    end
        
    if ~isempty(dir(fullfile(basepath,[basename,'.lfp']))) && ...
            ~isempty(dir(fullfile(basepath,'anatomical_map.csv')))
        
        % updates basename.session with channel map
         channel_mapping('basepath',basepath,'fig',false)
         
        session = loadSession(basepath,basename);
        try
            bad_channels = session.channelTags.Bad.channels;
        catch
            warning('no bad channels found. Consider verifying.')
            bad_channels = [];
        end
        
        channel = [];
        % Find cortical channel
        try
            cortex_channels = session.brainRegions.Cortex.channels;
            if ~isempty(bad_channels)
                cortex_channels(ismember(cortex_channels,bad_channels)) = [];
            end
            channel = cortex_channels(1);
            
         catch
            disp('Cortex not found, skipping session')
            continue
        end
        
          % use Dentate channel instead
        if isempty('channel')
            disp('No Cortex channel. Using last dentate channel')
            dentate_channels = session.brainRegions.Dentate.channels;
            dentate_channels(ismember(dentate_channels,bad_channels)) = [];
            channel = dentate_channels(end);
        end
        
            interictal_spikes = FindRipples('basepath',basepath,...
                'channel',channel,...
                'thresholds',[.25 3],...
                'durations',[10 150],...
                'minDuration',10,...
                 'passband',[40 250],...
                 'saveMat',false,...
                 'EMGThresh',.90); % Removes events that are above 90% corrrelated across all channels
             
            save(fullfile(basepath,[basename '.interictal_spikes.events.mat']),'interictal_spikes')
        
%         % Load ripples, if ripple event overlaps with IED then define as
%         % IED ( Gelinas et al 2016)
%         load(fullfile(basepath,[basename,'.ripples.events.mat']))
%                 
%         % Find IEDs that were also detected as physiological ripples
%        [status,interval,~] = InIntervals(interictal_spikes.peaks,ripples.timestamps);
%         
%        % Set those phyviological ripples as flagged 
%        ripples.flagged = interval(status);
%        
%        save(fullfile(basepath,[basename '.ripples.events.mat']),'ripples')

    end
end
end