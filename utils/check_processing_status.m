function check_processing_status(data_folder,varargin)
% Check df of basepaths for processing steps required for analysis 
%
% input: 
%   df: csv containing basepaths for checking

% output: 
%   saves session_status.csv to savepath indicating missing items from
%   basepath with 1 being found and 0 being not found. 

% Items that are checked: 
% Lfp file
% Sleep states
% Tracking: DLC output + animal behavior file + restricxy file + maze
% coords, restrictxy
% Kilosort folder (with rez file)
% Cell metrics
% restricxy file
% ripples
% ol_scoring
warning('off','all')

p = inputParser;
p.addParameter('check_lfp',true,@islogical)
p.addParameter('check_sleep',true,@islogical)
p.addParameter('check_tracking',true,@islogical)
p.addParameter('check_sorting',true,@islogical)
p.addParameter('check_cell_metrics',true,@islogical)
p.addParameter('check_ripples',true,@islogical)
p.addParameter('check_anatomical',true,@islogical)
p.addParameter('check_scoring',true,@islogical)
p.addParameter('remedy',true,@islogical)

p.parse(varargin{:})
check_lfp = p.Results.check_lfp;
check_sleep = p.Results.check_sleep;
check_tracking = p.Results.check_tracking;
check_sorting = p.Results.check_sorting;
check_cell_metrics = p.Results.check_cell_metrics;
check_ripples = p.Results.check_ripples;
check_anatomical = p.Results.check_anatomical;
check_scoring = p.Results.check_scoring;
remedy = p.Results.remedy;

% initialize session_status.csv as table
col_names = {'basepath','lfp','sleep_states','tracking_dlc',...
    'tracking_animalBehavior','tracking_restrictxy','tracking_mazeCoords',...
    'sorting_Kilosort','sorting_phyRez','tracking_sessionBehavioralTracking',...
    'cell_metrics','ripples','anatomical_map','behavior_scoring'};
session_status = cell2table(cell(0,length(col_names)),'VariableNames', col_names);

% find all sessions in data_folder
% handle input as csv of basepaths or directory
if contains(data_folder,'.csv')
    df = readtable(data_folder);
    data_folder = fileparts(data_folder);
else
    df = compile_sessions(data_folder);
end

basepaths = unique(df.basepath);
% loop through basepaths
for i = 1:length(basepaths)
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    
    session_status.basepath{i} = basepath;
    folder_contents = dir(basepath);
    
    if check_lfp
        session_status.lfp(i) = isfile(fullfile(basepath,[basename,'.lfp']));
    end
    
    if check_sleep
        session_status.sleep_states(i) = isfile(fullfile(basepath,[basename,'.SleepState.states.mat'])) ...
        & isfile(fullfile(basepath,[basename,'.SleepScoreLFP.LFP.mat'])) ...
        & isfile(fullfile(basepath,[basename,'.SleepStateEpisodes.states.mat']));
    end
    
    if check_tracking
        session_status.tracking_animalBehavior(i) = isfile(fullfile(basepath,[basename,'.animal.behavior.mat']));
        session_status.tracking_restrictxy(i) = isfile(fullfile(basepath,[basename,'.restrictxy.mat']));
        session_status.tracking_mazeCoords(i) = sum(isfile(fullfile(basepath,'*maze_coords.csv'))) >= 1;
        if isfile(fullfile(basepath,[basename,'.session.mat']))
            session = loadSession(basepath,basename);
            session_status.tracking_sessionBehavioralTracking(i) = isfield(session,'behavioralTracking');
        else
            session_status.tracking_sessionBehavioralTracking(i) = false;
        end
        vid_files = dir(fullfile(basepath,['*','.avi']));
        if length(find(contains({folder_contents.name},'DLC')))/length(vid_files) >= 3 % there are min 3 dlc files per video
            session_status.tracking_dlc(i) = true;
        else
            session_status.tracking_dlc(i) = false;
        end
    end
    
    if check_sorting
        
        if ~isempty(find(contains({folder_contents.name},'Kilosort')))
            ks = folder_contents(find(contains({folder_contents.name},'Kilosort'))).name;
            session_status.sorting_Kilosort(i) = true;
            % check if phy log has been created, thereby the user opened
            % session in phy
            session_status.sorting_phyRez(i) = isfile(fullfile(basepath,[ks,filesep,'phy.log']));
        else
            session_status.sorting_Kilosort(i) = false;
            session_status.sorting_phyRez(i) = false;
        end
    end
    
    if check_cell_metrics
        session_status.cell_metrics(i) = isfile(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']));
    end
    
    if check_ripples
        session_status.ripples(i) = isfile(fullfile(basepath,[basename,'.ripples.events.mat']));
    end
    
    if check_anatomical
        session_status.anatomical_map(i) = isfile(fullfile(basepath,'anatomical_map.csv'));
    end
    
    if check_scoring
        % if they have a video
        if sum(contains({folder_contents.name},'.avi')) >= 1
            session_status.behavior_scoring(i) = isfile(fullfile(basepath,[basename,'_ol_scoring.csv']));
        else
            session_status.behavior_scoring(i) = nan;
        end
    end
end

% save table to data_folder as csv
writetable(session_status,fullfile(data_folder,['session_check_',datestr(date),'.csv']))


end



