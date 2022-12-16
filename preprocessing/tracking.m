classdef tracking
    % functions for processing tracking, either from raw DLC output or
    % general behavior file
    
    methods(Static)
        
        function vec_out = scale_transform(vec,origin,scale_factor)
            % translate coordiante vector relative to origin and scale
            vec_out =(vec - origin)/scale_factor;
        end
        
        function ep_idx = get_epoch_idx(session,behavior)
            % for all behavior epochs, creates logical idx corresponding to epoch
            % boundaries. arrays are nested in cell corresponding to each behavioral
            % epoch (saved in same position as session.epochs)
            %
            % loop through epochs to retrieve start/end used in restrictxy below
            for ep = 1:length(session.behavioralTracking)
                epoch = session.behavioralTracking{1,ep}.epoch;
                
                ep_idx{ep} = (behavior.timestamps >= session.epochs{1, epoch}.startTime)...
                    & (behavior.timestamps <= session.epochs{1, epoch}.stopTime);
            end
            
        end
        
        function [x_origin,y_origin,scale_factor_x,scale_factor_y] = get_scale_params(session,behavior,basepath)
            % get origin and scale_factor from maze_coords.csv for each epoch. Depends on
            % maze_size being populated in session.behavioralTracking.maze_size (in cm)
            %
            % outputs are cell arrays for each behavioral epoch
            
            
            % loop through epochs to retrieve start/end used in restrictxy below
            for ep = 1:length(session.behavioralTracking)
                epoch = session.behavioralTracking{1,ep}.epoch;
                
                % load maze coords for given epoch video
                maze_coords_df = readtable(fullfile(basepath,...
                    [extractBefore(session.behavioralTracking{1,ep}.notes,'.avi'),'_maze_coords.csv']));
                
                % get max/min
                x_max = max(maze_coords_df.x(ismember(maze_coords_df.object,'corner')));
                x_min = min(maze_coords_df.x(ismember(maze_coords_df.object,'corner')));
                y_max = max(maze_coords_df.y(ismember(maze_coords_df.object,'corner')));
                y_min = min(maze_coords_df.y(ismember(maze_coords_df.object,'corner')));
                
                
                % compute parameters required for scaling
                x_origin{ep} = median(x_min:x_max);
                y_origin{ep} = median(y_min:y_max);
                
                scale_factor_x{ep} = (x_max - x_min)/behavior.epochs{1, epoch}.maze_size; %pixels/cm
                scale_factor_y{ep} = (y_max - y_min)/behavior.epochs{1, epoch}.maze_size; %pixels/cm
                
            end
            
        end
        
        function behavior = scale_coords(session,behavior,basepath)
            % updates behavior file scale and translate coordinates (centered at 0,0)
            
            
            % get max/min for maze per epoch from maze_coords.csv. Uses info from
            % session.behaviorTracking to load relevant maze_coords for each epoch.
            [x_origin,y_origin,scale_factor_x,scale_factor_y] = tracking.get_scale_params(session,behavior,basepath);
            idx = tracking.get_epoch_idx(session,behavior);
            
            % rescale coordinates relative to parameters for each epoch
            coord_names = fieldnames(behavior.position);
            
            for ep = 1:length(session.behavioralTracking)
                idx_ep = idx{ep};
                
                for i = find(contains(coord_names,{'x'}))'
                    x = behavior.position.(coord_names{i});
                    behavior.position.(coord_names{i})(idx_ep) = tracking.scale_transform(x(idx_ep),...
                        x_origin{ep},scale_factor_x{ep});
                end
                
                for i = find(contains(coord_names,{'y'}))'
                    y = behavior.position.(coord_names{i});
                    behavior.position.(coord_names{i})(idx_ep) = tracking.scale_transform(y(idx_ep),...
                        y_origin{ep},scale_factor_y{ep});
                end
                
                % update speed to match new units
                [speed_ep, accel_ep,~] = linear_motion(behavior.position.x(idx_ep),...
                    behavior.position.y(idx_ep),behavior.sr,.1);
                behavior.speed(idx_ep) = speed_ep;
                behavior.acceleration(idx_ep) = accel_ep;
            end
            
            behavior.position.units = 'cm';
            
        end
        
        function behavior = restrict(session,behavior,basepath)
            % input:
            %   session: CellExplorer session file with completed epochs,
            %       and behavioralTracking.
            %   behavior: genearl behavior file
            % output:
            %   updated behavior file with position data restricted
            %   relative to .restrictxy.mat file (created with
            %   manual_trackerjumps)
            basename = basenameFromBasepath(basepath);
            
            % get start and stop from session epochs to use for limited
            % coordinates to given epochs
            start = [];
            stop = [];
            % gets index to restrict xy
            for ep = 1:length(session.behavioralTracking)
                epoch = session.behavioralTracking{1,ep}.epoch;
                start = [start,session.epochs{epoch}.startTime];
                stop = [stop,session.epochs{epoch}.stopTime];
                
            end
            
            % Get good index (coordinates within boundary) from file or
            % made with manual_trackerjumps.
            if ~isempty(dir(fullfile(basepath,[basename,'.restrictxy.mat'])))
                load(fullfile(basepath,[basename,'.restrictxy.mat']))
            else
                good_idx = manual_trackerjumps(behavior.timestamps,...
                    behavior.position.x,...
                    behavior.position.y,...
                    start,...
                    stop,...
                    basepath,'darkmode',false);
            end
            % for primary xy coordinates, restrict to boundaries provided by good_idx
            behavior.position.x(~good_idx) = NaN;
            behavior.position.y(~good_idx) = NaN;
        end
        
        function [t,x,y,v,a,angle,units,source,fs,notes,extra_points,vidnames] = ...
                extract_tracking(basepath,primary_coords_dlc,likelihood_dlc,smooth_factor)
            
            % initalize variables to pull
            t = [];
            x = [];
            y = [];
            v = [];
            a = [];
            angle = [];
            units = [];
            source = [];
            notes = [];
            extra_points = [];
            vidnames = [];
            
            % below are many methods on locating tracking data from many formats
            [tracking_struct,field_names] = tracking.process_and_sync_dlc_SNLab('basepath',basepath,...
                'primary_coords',primary_coords_dlc,...
                'likelihood',likelihood_dlc);
            
            t = tracking_struct.timestamps;
            fs = 1/mode(diff(t));
            vidnames = tracking_struct.vidnames;
            
            x = tracking_struct.position.x(:,primary_coords_dlc);
            y = tracking_struct.position.y(:,primary_coords_dlc);
            
            if length(primary_coords_dlc) > 1
                angle = xy_angle(x,y);
                % compute average point between two coords
                x = nanmedian(x,2);
                y = nanmedian(y,2);
                [v, a,~] = linear_motion(x,y,fs,smooth_factor);
            else
                [v, a,~] = linear_motion(x,y,fs,smooth_factor);
            end
            
            % multiple tracking points will likely exist, extract here
            x_col = field_names(contains(field_names,'x'));
            y_col = field_names(contains(field_names,'y'));
            extra_points = struct();
            for i = 1:length(x_col)
                extra_points.([x_col{i},'_point']) = tracking_struct.position.x(:,i);
                extra_points.([y_col{i},'_point']) = tracking_struct.position.y(:,i);
            end
            
            units = 'pixels';
            source = 'deeplabcut';
            
            if length(t) > length(x)
                t = t(1:length(x));
            elseif length(x) > length(t)
                x = x(1:length(t));
                y = y(1:length(t));
                % adjust other tracker points
                for name = fields(extra_points)'
                    extra_points.(name{1}) = extra_points.(name{1})(1:length(t));
                end
            end
            
            notes = ['primary_coords: ',num2str(primary_coords_dlc),...
                ', likelihood: ',num2str(likelihood_dlc)];
            
        end
        
        % Sync dlc to ttl for ephys
        function [tracking_struct,field_names] = process_and_sync_dlc_SNLab(varargin)
            % Unpacks DLC and syncs with digitalin.events.mat timestamps
            %
            % Run this after you have exported deeplabcut csv results, have
            % basename.session.epochs labeled, and have digitalin.events.mat.
            %
            % optional inputs:
            % basepath: location of SNLab session data (contains basename.session file,
            %  digitalin.events.mat, and DLC output saved as CSV).
            % video_channel_ttl: digitalin channel that contains video ttl pulses (default channel 1).
            % primary_coords: cordinates of interest (default first dlc bodypart)
            %
            %
            % L Berkowitz 03/2022 adapted some code by R Harvey
            
            % to-do:
            % add functionality to choose primary coordinates
            
            % parse inputs
            p = inputParser;
            p.addParameter('basepath',pwd,@isfolder);
            p.addParameter('video_channel_ttl',1,@isnumeric); % digitalin channel that contains video ttl
            p.addParameter('primary_coords',1:2,@isnumeric); % which tracker point do you want
            p.addParameter('likelihood',.95,@isnumeric); % tracking quality thres [0-1]
            p.addParameter('pulses_delta_range',0.01,@isnumeric); % range for ttls
            
            p.parse(varargin{:});
            primary_coords = p.Results.primary_coords; % not used currently
            basepath = p.Results.basepath;
            video_channel_ttl = p.Results.video_channel_ttl;
            likelihood = p.Results.likelihood;
            pulses_delta_range = p.Results.pulses_delta_range;
            basename = basenameFromBasepath(basepath);
            
            % check for events, cannot process without them
            if exist(fullfile(basepath,'digitalin.events.mat'),'file') % only load if timestamps have been processed
                load(fullfile(basepath,'digitalin.events.mat'));
                
                if exist('parsed_digitalIn','var')
                    digitalIn = parsed_digitalIn;
                end
            else % if none, make them and update basename.session
                disp('processing events, one moment')
                % make digitalin.events.mat
                process_digitalin(basepath,session.extracellular.sr)
                % update epochs
                update_epochs('basepath',basepath,'annotate',true)
                % update behavioralTracking
                update_behavioralTracking('basepath',basepath)
            end
            
            % check for behavioralTracking field from session file and if present, grab
            % tracking info (tracking file name, epoch index, and frame
            % rate
            
            % first load session file
            session = loadSession(basepath,basename);
            
            if isfield(session,'behavioralTracking')
                
                [dlc_files,vidnames,ep_idx,frame_rate] =  tracking.get_tracking_info_from_session(session,basename);
                
            else % create it and reload session and pull dlc,video, and epoch info
                % update session with dlc/video info
                update_behavioralTracking('basepath',basepath)
                session = loadSession(basepath,basename); % reload updated session file
                % pull tracking info
                [dlc_files,vidnames,ep_idx,frame_rate] =  tracking.get_tracking_info_from_session(session,basename);
                
            end
            
            % grab video ttls
            behav = 1;
            for epoch = ep_idx
                % grab the video ttl timestamps within the epoch boundaries
                idx = digitalIn.timestampsOn{1, video_channel_ttl} >= session.epochs{1, epoch}.startTime  & ...
                    digitalIn.timestampsOn{1, video_channel_ttl} <= session.epochs{1, epoch}.stopTime;
                video_ttl{behav} = digitalIn.timestampsOn{1, video_channel_ttl}(idx);
                behav = behav + 1;
                disp(['Epoch ', num2str(epoch), ' contains behavior tag. Grabbing video channel ttls for this epoch.'])
            end
            
            %% Sync video ttl with trackin file
            for ii = 1:length(dlc_files)
                disp(['processing file ',num2str(ii), ' of ',num2str(length(dlc_files))])
                
                % get frame rate of video
                fs = frame_rate(ii);
                
                % load csv with proper header
                opts = detectImportOptions(fullfile(basepath,dlc_files{ii}),'NumHeaderLines',2);
                df = readtable(fullfile(basepath,dlc_files{ii}),opts);
                
                % get names of fields, these will be as long as tracker points
                % used times 3 because [x,y,likelihood]
                field_names = fields(df);
                
                % locate columns with [x,y,likelihood]
                x_col = find(contains(field_names,'x'));
                y_col = find(contains(field_names,'y'));
                likelihood_col = find(contains(field_names,'likelihood'));
                
                % filter out bad tracker points by likelihood thres
                for i = 1:length(x_col)
                    idx = df{:,likelihood_col(i)} < likelihood;
                    df{idx,x_col(i)} = NaN;
                    df{idx,y_col(i)} = NaN;
                end
                ts = df{:,1}/fs;
                x = df{:,x_col};
                y = df{:,y_col};
                
                % store tracking for each video file
                tempTracking{ii} = tracking.sync_ttl(basepath,video_ttl{ii},x,y,ts,fs,pulses_delta_range);
                trackFolder(ii) = ii;
            end
            
            % Concatenating tracking fields...
            x = []; y = []; timestamps = []; folder = []; samplingRate = []; description = [];
            for ii = 1:size(tempTracking,2)
                x = [x; tempTracking{ii}.position.x];
                y = [y; tempTracking{ii}.position.y];
                timestamps = [timestamps; tempTracking{ii}.timestamps];
                folder{ii} = tempTracking{ii}.folder;
                samplingRate = [samplingRate; tempTracking{ii}.samplingRate];
                description{ii} = tempTracking{ii}.description;
            end
            
            % save data to ouput structure
            tracking_struct.position.x = x;
            tracking_struct.position.y = y;
            tracking_struct.folders = folder;
            tracking_struct.samplingRate = samplingRate;
            tracking_struct.timestamps = timestamps;
            tracking_struct.description = description;
            tracking_struct.vidnames = vidnames;
            
        end
       
        function [tracking_struct] = sync_ttl(basepath,video_ttl,x,y,ts,fs,pulses_delta_range)
            
            % if ~exist(fullfile(folder,'digitalIn.events.mat'),'file')
            %     digitalIn = getDigitalIn('all','folder',folder);
            % end
            %
            % load(fullfile(folder,'digitalIn.events.mat'))
            %
            % Len = cellfun(@length, digitalIn.timestampsOn, 'UniformOutput', false);
            % [~,idx] = max(cell2mat(Len));
            % bazlerTtl = digitalIn.timestampsOn{idx};
            
            %check for extra pulses of much shorter distance than they should
            extra_pulses = diff(video_ttl)<((1/fs)-(1/fs)*pulses_delta_range);
            video_ttl(extra_pulses) = [];
            
            video_intan_diff = length(video_ttl) - size(x,1);
            
            [x,y,ts,video_ttl] = tracking.match_video_frames_to_ttl(video_ttl,video_intan_diff,x,y,ts,fs);
            
            [~,folder_name] = fileparts(basepath);
            tracking_struct.position.x = x;
            tracking_struct.position.y = y;
            tracking_struct.timestamps = video_ttl;
            tracking_struct.originalTimestamps = ts;
            tracking_struct.folder = folder_name;
            tracking_struct.samplingRate = fs;
            tracking_struct.description = '';
        end
        
        function [x,y,t,video_ttl] = match_video_frames_to_ttl(video_ttl,basler_intan_diff,x,y,t,fs)
            
            % match basler frames con ttl pulses
            if (length(video_ttl) == size(x,1)) || abs(basler_intan_diff)<=2 %assumes 1 frame could be cut at 0 and 1 frame at end
                disp('N of frames match!!');
            elseif basler_intan_diff>0 && abs(basler_intan_diff)<fs
                disp([num2str(abs(length(video_ttl) - size(x,1))) ' of frames dont match, probably at the end of the recording']);
                video_ttl = video_ttl(1:size(x,1));
            elseif basler_intan_diff<0 && abs(basler_intan_diff)<fs
                disp([num2str(abs(length(video_ttl) - size(x,1))) ' of frames dont match, probably at the beggining of the recording']);
                x = x(1:length(video_ttl),:);
                y = y(1:length(video_ttl),:);
            elseif basler_intan_diff<0 && abs(basler_intan_diff)>fs
                disp([num2str(abs(size(x,1) - length(video_ttl)))...
                    ' video frames without TTL... was the recording switched off before the camera? Cutting positions accordingly...']);
                x = x(1:length(video_ttl),:);
                y = y(1:length(video_ttl),:);
            elseif abs(basler_intan_diff)>2*fs
                warning('More than 2 seconds missalignment in total in this session...will adjust to the closer one...');
                if basler_intan_diff>0
                    video_ttl = video_ttl(1:size(x,1));
                else
                    x = x(1:length(video_ttl),:);
                    y = y(1:length(video_ttl),:);
                end
            elseif isempty(video_ttl)
                video_ttl = t;
            else
                warning('Unnoticed problem with Camera/Intan... I would go back and check both step by step');
            end
        end
        
        function [tracking_files,vid_names,ep_idx,frame_rate] =  get_tracking_info_from_session(session,basename)
            % extracts tracking information (tracking filename, video name, and epoch index
            % basename.session.behavioralTracking
            
            disp(['Checking for DLC files in ', basename, '.session.behavioralTracking'])
            
            for i = 1:length(session.behavioralTracking)
                tracking_files{i} = session.behavioralTracking{1,i}.filenames;
                vid_names{i} = session.behavioralTracking{1,i}.notes;
                ep_idx(i) =  session.behavioralTracking{1,i}.epoch;
                frame_rate(i) = session.behavioralTracking{1,i}.framerate;
            end
            
        end
        
        % General  
        function mdl_params = grab_dlc_crop(dlc_file_name,varargin)
            % pulls cropping parameters from dlc_crop_path (csv with croppsing
            % parameters for model indicated in dlc output.
            p=inputParser;
            p.addParameter('dlc_crop','C:\Users\schafferlab\github\SNLab_ephys\behavior\dlc_crop_parameters.csv',@ischar);
            parse(p,varargin{:});
            
            dlc_crop_path = p.Results.dlc_crop;
            
            % load csv with model parameters for cropping
            dlc_crop = readtable(dlc_crop_path);
            
            % initialize output
            mdl_params = table;
            % loop through models
            for i = 1:length(dlc_crop.model)
                % save parameters if file name contains the model name
                if contains(dlc_file_name,dlc_crop.model{i})
                    mdl_params = dlc_crop(i,:);
                end
            end
            
        end
        
    end
end