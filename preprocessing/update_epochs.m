function update_epochs(varargin)
% update_epochs will add epochs to basename.session file using
% digitalin.events.mat file.

% input:
%   basepath: path to session folder
%   acquisition_event_flag: True if digitalin.events contain start and end
%   epcohs for whole session (i.e. start, start, end, end for start
%   session, start video, end video, end session). This approach was
%   discontinued in SNLab but appears in some data.
%   epoch_events: index for digitalIn.events.mat for epoch related ttls
%   annotate: Option for adding details to epochs via gui_session

% output:
%   basename.session with updated epochs

% input parser
% Outputs
p=inputParser;
addParameter(p,'basepath',pwd); % single or many basepaths in cell array or uses pwd
addParameter(p,'acquisition_event_flag',false); % overwrite previously saved data (will remove custom fields)
addParameter(p,'epoch_events',2); % index for where epoch events are saved in digitalIn.events.mat file
addParameter(p,'overwrite',true); % by default update epochs
addParameter(p,'annotate',false); % save animal.behavior.mat

% addParameter(p,'maze_size',30); % maze size in cm

parse(p,varargin{:});
basepath = p.Results.basepath;
acquisition_event_flag = p.Results.acquisition_event_flag;
annotate = p.Results.annotate;
overwrite = p.Results.overwrite;
epoch_events = p.Results.epoch_events;
basename = basenameFromBasepath(basepath);

if exist(fullfile(basepath,'digitalIn.events.mat'),'file')
    
    % load digitalin file
    load(fullfile(basepath,'digitalIn.events.mat'),'digitalIn')
    % load session file
    load(fullfile(basepath,[basename, '.session.mat']),'session')
    
    % check if session.epochs were already created
    if rem(length(digitalIn.timestampsOn{1, epoch_events}),length(session.epochs)) == 1
        % if length of events are not divisible by n epochs, visuallly
        % inspect epochs
        disp('N epochs does not match number of timestamps. Check epochs and indicate if they''ve been added')
        gui_session(basepath)
        overwrite = ui_check_epochs;
    end
    
    % if epochs haven't been added, add them here
    if overwrite
        if acquisition_event_flag
            start_idx = 2;
            % first and last time stamp are always acquisition
            session.epochs{1}.name =  'acquisition';
            session.epochs{1}.startTime =  parsed_digitalIn.timestampsOn{1, 2}(1);
            session.epochs{1}.stopTime = parsed_digitalIn.timestampsOn{1, 2}(end);
            
        else
            start_idx = 1;
        end
        
        % loop through the other epochs
        ii = start_idx;
        for i = start_idx:2:size(parsed_digitalIn.timestampsOn{1, 2},1)-1 % by default 2nd column is events
            session.epochs{ii}.name =  char(i);
            session.epochs{ii}.startTime =  parsed_digitalIn.timestampsOn{1, 2}(i);
            session.epochs{ii}.stopTime =  parsed_digitalIn.timestampsOff{1, 2}(i+1);
            ii = ii+1;
        end
    end
    
    % save updated session file
    save(fullfile(basepath,[basename, '.session.mat']),'session');
    
    % if user wants to annotate epoch details, open the gui
    if annotate
        % verify videos added
        disp('check epoch notes and verify videos were added to notes')
        gui_session(basepath)
    end
else
    warning('No digitalin events found in basepath. Epochs not updated in session file.')
    return
end
end

function run_update_epoch = ui_check_epochs
%
answer = questdlg('Have epochs been added?', ...
    'Update epochs', ...
    'yes','no','no');

% Handle response
switch answer
    case 'yes'
        disp([answer ' epochs have already been added to session file'])
        run_update_epoch = false;
    case 'no'
        disp([answer ' epochs not been added to session file. Adding them now'])
        run_update_epoch = true;
end

end