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

p.parse(varargin{:})
check_lfp = p.Results.check_lfp;
check_sleep = p.Results.check_sleep;
check_tracking = p.Results.check_tracking;
check_sorting = p.Results.check_sorting;
check_cell_metrics = p.Results.check_cell_metrics;
check_ripples = p.Results.check_ripples;
check_anatomical = p.Results.check_anatomical;
check_scoring = p.Results.check_scoring;

% find all sessions in data_folder
% handle input as csv of basepaths or directory
if contains(data_folder,'.csv')
    df = readtable(data_folder);
    data_folder = fileparts(data_folder);
    basepaths = unique(df.basepath);
else
    df = compile_sessions(data_folder);
    basepaths = unique([df.basepath{:}]);
end

% initialize table 
% initialize session_status.csv as table
col_names = {'basepath','lfp','lfp_date',...
    'sleep_states','sleep_states_date',...
    'tracking_dlc','tracking_dlc_date',...
    'tracking_animalBehavior','tracking_animalBehavior_date',...
    'tracking_restrictxy','tracking_restrictxy_date',...
    'sorting_Kilosort', 'sorting_Kilosort_date',...
    'sorting_phyRez','sorting_phyRez_date',...
    'tracking_sessionBehavioralTracking','tracking_sessionBehavioralTracking_date',...
    'cell_metrics','cell_metrics_date',...
    'ripples','ripples_date',...
    'anatomical_map','anatomical_map_date',...
    'behavior_scoring','behavior_scoring_date'};

col_types = {'string','logical','string',...
    'logical','string',...
    'logical','string',...
    'logical','string',...
    'logical','string',...
    'logical', 'string',...
    'logical','string',...
    'logical','string',...
    'logical','string',...
    'logical','string',...
    'logical','string',...
    'logical','string'};

session_status = table('Size',[length(basepaths), length(col_names)],'VariableNames',col_names,'VariableTypes',col_types);

% loop through basepaths

for i = 1:length(basepaths)
    basepath = basepaths{i};
    basename = basenameFromBasepath(basepath);
    
    session_status.basepath{i} = basepath;
    folder_contents = dir(basepath);
    
    if check_lfp
        nameExt=fullfile(basepath,[basename,'.lfp']);
        session_status.lfp(i) = isfile(nameExt);

        if  ismember(session_status.lfp(i),1)
            %session_status.lfp_date(i)=dir(fullfile(basepath,[basename,'.lfp'])).date;
            session_status.lfp_date(i)=insert_date(nameExt);
        end

    end
    
    if check_sleep
        nameExt=fullfile(basepath,[basename,'.SleepState.states.mat']);
        session_status.sleep_states(i) = isfile(nameExt) ...
        & isfile(fullfile(basepath,[basename,'.SleepScoreLFP.LFP.mat'])) ...
        & isfile(fullfile(basepath,[basename,'.SleepStateEpisodes.states.mat']));
        
        if  ismember(session_status.lfp(i),1)
            session_status.sleep_states_date(i)=insert_date(nameExt);
        end

    end
    
    if check_tracking
        nameExt=fullfile(basepath,[basename,'.animal.behavior.mat']);
        session_status.tracking_animalBehavior(i) = isfile(nameExt);
        if  ismember(session_status.tracking_animalBehavior(i),1)
            session_status.tracking_animalBehavior_date(i)=insert_date(nameExt);
        end
        nameExt=fullfile(basepath,[basename,'.restrictxy.mat']);
        session_status.tracking_restrictxy(i) = isfile(nameExt);
        if  ismember(session_status.tracking_animalBehavior(i),1)
            session_status.tracking_restrictxy_date(i)=insert_date(nameExt);
        end
        if isfile(fullfile(basepath,[basename,'.session.mat']))
            session = loadSession(basepath,basename);
            session_status.tracking_sessionBehavioralTracking(i) = isfield(session,'behavioralTracking');
            if ismember(session_status.tracking_sessionBehavioralTracking(i),1)
                sessionName=fullfile(session.general.basePath, ...
                    session.behavioralTracking{1,1}.filenames);
                session_status.tracking_sessionBehavioralTracking_date(i) = insert_date(sessionName);
            end
        else
            session_status.tracking_sessionBehavioralTracking(i) = false;
        end
        
        vid_files =fullfile(basepath,['*','.avi']);
        if length(find(contains({folder_contents.name},'DLC')))/length(dir(vid_files)) >= 3 % there are min 3 dlc files per video
            session_status.tracking_dlc(i) = true;
            if  ismember(session_status.tracking_dlc(i),1)
                
                
                session_status.tracking_dlc_date(i)=date_compare(folder_contents,'DLC');
                %TODO What if theres multiple videos
            end
        else
            session_status.tracking_dlc(i) = false;
        end
    end
    
    if check_sorting
        kiloPresent=find(contains({folder_contents.name},'Kilosort'));
        if ~isempty(kiloPresent)
            ks = folder_contents(find(contains({folder_contents.name},'Kilosort'))).name;
            session_status.sorting_Kilosort(i) = true;
            %TODO when KiloPresent has two kilofolders
            session_status.sorting_Kilosort_date(i)=date_compare(folder_contents,'Kilosort');
            nameExt=fullfile(basepath,[ks,filesep,'phy.log']);
            % check if phy log has been created, thereby the user opened
            % session in phy
            
            session_status.sorting_phyRez(i) = isfile(nameExt);
            if ismember(session_status.sorting_phyRez(i),1)
                session_status.sorting_phyRez_date(i)=insert_date(nameExt);
            end
   
        else
            session_status.sorting_Kilosort(i) = false;
            session_status.sorting_phyRez(i) = false;
        end

    end
    
    if check_cell_metrics
        nameExt=fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']);
        session_status.cell_metrics(i) = isfile(nameExt);
        if ismember(session_status.cell_metrics(i),1)
            session_status.cell_metrics_date(i)=insert_date(nameExt);
        end
    end
    
    if check_ripples
        nameExt=fullfile(basepath,[basename,'.ripples.events.mat']);
        session_status.ripples(i) = isfile(nameExt);
        if ismember(session_status.ripples(i),1)
            session_status.ripples_date(i)=insert_date(nameExt);
        end
    end
    
    if check_anatomical
        nameExt=fullfile(basepath,'anatomical_map.csv');
        session_status.anatomical_map(i) = isfile(nameExt);
        if ismember(session_status.anatomical_map(i),1)
            session_status.anatomical_map_date(i)=insert_date(nameExt);
        end
    end
    
    if check_scoring
        % if they have a video
        if sum(contains({folder_contents.name},'.avi')) >= 1
            nameExt=fullfile(basepath,[basename,'_ol_scoring.csv']);    
            session_status.behavior_scoring(i) = isfile(nameExt);
            if ismember(session_status.behavior_scoring(i),1)
                session_status.behavior_scoring_date(i)=insert_date(nameExt);
            end
        
        else
            session_status.behavior_scoring(i) = 0;  %TODO CHECK
        end
    end
end

% save table to data_folder as csv
writetable(session_status,fullfile(data_folder,['session_check_',datestr(date),'.csv']))


end
%Get date into proper format
function out=date_compare(folder_contents,key)
    names=folder_contents(find(contains({folder_contents.name},key)));
    out=datetime(names(1).date);
    if length(names)>1
        for i=2:length(names)
            if datetime(names(i).date)>out
            out=datetime(names(i).date);
            end
        end
    end
    out=datestr(out,'mm/dd/yy');
end
function out=insert_date (name)
    tempPath=dir(name);
    out=datestr(tempPath(length(tempPath)).date, 'mm/dd/yy');
end




