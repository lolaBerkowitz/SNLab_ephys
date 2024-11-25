classdef tracking
    % functions for processing tracking, either from raw DLC output or
    % general behavior file
    
    methods(Static)
        
       
        function vec_out = scale_transform(vec,origin,scale_factor)
            % translate coordiante vector relative to origin and scale
            vec_out =(vec - origin)/scale_factor;
        end
        
        function vec_out = scale(vec,scale_factor)
            
            vec_out = vec/scale_factor;
            
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
            
            basename = basenameFromBasepath(basepath);
            
            % loop through epochs to retrieve start/end used in restrictxy below
            for ep = 1:length(session.behavioralTracking)
                
                % exclude VR 
                if contains(session.behavioralTracking{1,ep}.filenames,'godot')
                    continue
                end
                
                epoch = session.behavioralTracking{1,ep}.epoch;
                
                pixel_distance = maze_distance_gui(video_path,maze_sizes);
                
                % load maze coords for given epoch video
                maze_coords_df = readtable(fullfile(basepath,...
                    [extractBefore(session.behavioralTracking{1,ep}.notes,'.avi'),'_maze_coords.csv']));
                pixel_distance
                % get max/min
                x_max = max(maze_coords_df.x(ismember(maze_coords_df.object,'corner')));
                x_min = min(maze_coords_df.x(ismember(maze_coords_df.object,'corner')));
                y_max = max(maze_coords_df.y(ismember(maze_coords_df.object,'corner')));
                y_min = min(maze_coords_df.y(ismember(maze_coords_df.object,'corner')));
                
                
                % compute parameters required for scaling
                x_origin{ep} = median(x_min:x_max);
                y_origin{ep} = median(y_min:y_max);
                
                if length(behavior.epochs{1, epoch}.maze_size) > 1
                    scale_factor_x{ep} = (x_max - x_min)/behavior.epochs{1, epoch}.maze_size(1); %pixels/cm
                    scale_factor_y{ep} = (y_max - y_min)/behavior.epochs{1, epoch}.maze_size(2); %pixels/cm
                else 
                    scale_factor_x{ep} = (x_max - x_min)/behavior.epochs{1, epoch}.maze_size(1); %pixels/cm
                    scale_factor_y{ep} = (y_max - y_min)/behavior.epochs{1, epoch}.maze_size(1); %pixels/cm
                end
                % add scaled parameters to maze_coord_df
                maze_coords_df.x_scaled = tracking.scale_transform(maze_coords_df.x,x_origin{ep},scale_factor_x{ep});
                maze_coords_df.y_scaled = tracking.scale_transform(maze_coords_df.y,y_origin{ep},scale_factor_y{ep});
                
                % save scaled parameters to behavioralTracking 
                session.behavioralTracking{1, ep}.maze_coords = maze_coords_df;
           
                % save back to basepath
                writetable(maze_coords_df,fullfile(basepath,...
                    [extractBefore(session.behavioralTracking{1,ep}.notes,'.avi'),'_maze_coords.csv']))
            end
            save(fullfile(basepath,[basename,'.session.mat']),'session')
            
        end
        
        function behavior = scale_coords(session,behavior,basepath)
            basename = basenameFromBasepath(basepath);
            % updates behavior file scale and translate coordinates (centered at 0,0)
            if contains(behavior.position.units,'cm')
                disp('Coordinates already scaled')
                return
            end
            
            % get max/min for maze per epoch from maze_coords.csv. Uses info from
            % session.behaviorTracking to load relevant maze_coords for each epoch.
            [x_origin,y_origin,scale_factor_x,scale_factor_y] = tracking.get_scale_params(session,behavior,basepath);
            idx = tracking.get_epoch_idx(session,behavior);
            
            % rescale coordinates relative to parameters for each epoch
            coord_names = fieldnames(behavior.position);

            for ep = 1:length(session.behavioralTracking)
                
                % exclude VR 
                if contains(session.behavioralTracking{1,ep}.filenames,'godot')
                    continue
                end
                
                idx_ep = idx{ep};
                
                for i = find(contains(coord_names,{'x'}))'
                    x = behavior.position.(coord_names{i});
                    if length(idx_ep) ~= length(x)
                        continue
                    end
                    behavior.position.(coord_names{i})(idx_ep) = tracking.scale_transform(x(idx_ep),...
                        x_origin{ep},scale_factor_x{ep});
                end
                
                for i = find(contains(coord_names,{'y'}))'
                    y = behavior.position.(coord_names{i});
                    if length(idx_ep) ~= length(y)
                        continue
                    end
                    behavior.position.(coord_names{i})(idx_ep) = tracking.scale_transform(y(idx_ep),...
                        y_origin{ep},scale_factor_y{ep});
                end
               
                % BROKEN 
              % update speed to match new units 
              % TOO FIX - broken but not
              % necessarily needed as we just recreate them during analysis
              
%             [speed_, accel_,~] = linear_motion(behavior.position.x(idx_ep),...
%                     behavior.position.y(idx_ep),session.behavioralTracking{1, ep}.framerate,.1);
%             behavior.speed(idx_ep(2:end)) = speed_;
%             behavior.acceleration(idx_ep(2:end)) = accel_;
            end
            
            behavior.position.units = 'cm';
            save(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')

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
            
            save(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')

        end
        
        function [t,x,y,v,a,angle,units,source,fs,notes,extra_points,vidnames] = ...
                extract_godot_tracking(basepath)
            
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
            tracking_struct = tracking.process_and_sync_godot_SNLab('basepath',basepath);
            
            t = tracking_struct.timestamps;
            fs = nan;
            vidnames = tracking_struct.vidnames;
            v = tracking_struct.v; 
            
            x = tracking_struct.position.x*100;
            y = tracking_struct.position.y*100;
            angle = tracking_struct.angle; 
            
            % remove spurious velocities due to end of track
            v = smooth(v,round((1/median(diff(t)))/2)); 
            % calculate acceleration
            a = gradient(v);
       
            units = 'virtual_cm';
            source = 'godot';
            
            if length(t) > length(x)
                t = t(1:length(x));
            elseif length(x) > length(t)
                x = x(1:length(t));
                y = y(1:length(t));
            end
            
            notes = [];   
        end
        

        function pixel_distance = maze_distance_gui(video_path,maze_sizes,vid_time)
        % maze_distance_gui: manually get distance between points in pixels
        %
        % This fuction was made in order to get the ratio to convert tracking
        % points from pixels to cm.
        %
        % Input: 
        %   video_path: full path to video
        %
        % Output:
        %   pixel_distance: distance between clicked points in pixels
        %
        % Copyright (C) 2022 Ryan Harvey

        % check if video exists
        if ~exist(video_path,'file')
           error([video_path,'  video does not exist']) 
        end
        

%         % read video
%         vid_obj = VideoReader(video_path);
% 
%         % read first 10 seconds of video frames
%         frames = read(vid_obj, [1, round(vid_obj.FrameRate*1)]);
%         
        basepath = fileparts(video_path);
        img_path = fullfile(basepath,'temp_img.png');
        % create image save to current directory
        sys_cmd = ['ffmpeg -ss ', num2str(vid_time),' -i ',video_path,' -vframes 1 ',img_path];
        system(sys_cmd)
        
        im = imread(img_path); 
    
%         % init matrix to store flattened frames
%         grey_frames = zeros(vid_obj.Height, vid_obj.Width, size(frames,4));
%         for i = 1:size(frames, 4)
%             grey_frames(:,:,i) = rgb2gray(frames(:,:,:,i));
%         end
%         % take averge of 10 seconds
%         grey_frames_avg = mean(grey_frames, 3);

        % plot average frame
        fig = figure;
        imagesc(im) % display t he first frame        
        hold on;
        colormap('gray');
        axis('image')
        if nargin>1
            try
                title(['click on 2 key points approximately ' num2str(maze_sizes) 'cm apart and hit "enter" '])
            catch
                title(['click on 2 key points approximately ' num2str(maze_sizes) 'cm apart and hit "enter" '])
            end
        else
            title('click on 2 key points and hit "enter" ')
        end

        % let the user click around the coordinates
        corners = [];
        i = 1;
        while true
            [X,Y]=ginput(1);
            if isempty(X)
                break
            end
            corners = [corners;[X,Y]];
            plot(corners(:,1),corners(:,2),'r',X,Y,'*r')
            i=i+1;
        end
        close(fig)
        delete(img_path)

        % calculate pixel distance
        pixel_distance = pdist(corners,'euclidean');
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
        
        % Sync gadot to ttl for ephys 
        function [tracking_struct] = process_and_sync_godot_SNLab(varargin)
            
             % parse inputs
            p = inputParser;
            p.addParameter('basepath',pwd,@isfolder);
            p.addParameter('sync_ttl_channel',7,@isnumeric); % digitalin channel that contains video ttl           
            
            p.parse(varargin{:});
            basepath = p.Results.basepath;
            sync_ttl_channel = p.Results.sync_ttl_channel;
            
            basename = basenameFromBasepath(basepath);
            session = loadSession(basepath,basename);
            
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
            
            if isfield(session,'behavioralTracking')
                
                [godot_files,vidnames,ep_idx,~] =  tracking.get_tracking_info_from_session(session,basename);
                
            else % create it and reload session and pull dlc,video, and epoch info
                % update session with dlc/video info
                update_behavioralTracking('basepath',basepath)
                session = loadSession(basepath,basename); % reload updated session file
                % pull tracking info
                [godot_files,vidnames,ep_idx,~] =  tracking.get_tracking_info_from_session(session,basename);
                
            end
            ep_idx = ep_idx(contains(godot_files,'godot'));
            godot_files = godot_files(contains(godot_files,'godot'));
            for i = 1:length(godot_files)
                
                % load godot 
                vr_pos = readtable(fullfile(basepath,godot_files{i}));

               % get intan timestamp associated with first index 
                sync_ts = digitalIn.timestampsOn{1, sync_ttl_channel};
                % limit sync_ts to epoch
                keep_idx = sync_ts > session.epochs{1,ep_idx}.startTime & sync_ts < session.epochs{1,ep_idx}.stopTime;
                sync_ts = sync_ts(keep_idx);
                
                % find the first reward 
                reward_interval = findIntervals(logical(vr_pos.reward));
%                 between_reward_interval = findIntervals(~logical(vr_pos.reward));

                
                
                behav_ts = vr_pos.experiment_ts(reward_interval(1,1));
                computer_time_intan_diff = sync_ts(1) - behav_ts;
                vr_ts_sync = vr_pos.experiment_ts + computer_time_intan_diff;
                
                
%                 for interval_i = 1:size(reward_interval,1)
%                     
%                     behav_ts = vr_pos.experiment_ts(reward_interval(interval_i,1));
%                     computer_time_intan_diff(interval_i) = sync_ts(interval_i) - behav_ts;
%                     
%                 end
                
                
                % interpolate with the session time to get the remaining
                % timestamps 
%                 interp_times = fill_in_missing_timestamps(sync_ts,diff(vr_pos.experiment_ts)); 
                
                tempTracking{i}.position.x = vr_pos.z;
                tempTracking{i}.position.y = vr_pos.x;
                tempTracking{i}.v = sqrt(vr_pos.pitch.^2 + vr_pos.roll.^2);
                tempTracking{i}.angle = vr_pos.xz_angle;
                tempTracking{i}.timestamps = vr_ts_sync;
                tempTracking{i}.folder = vidnames{i};
                tempTracking{i}.samplingRate = NaN;
                tempTracking{i}.description = vidnames{i};
            end 
            
            % Concatenating tracking fields...
            x = []; y = []; v = []; xy_angle = []; timestamps = []; folder = []; samplingRate = []; description = [];
            for ii = 1:size(tempTracking,2)
                x = [x; tempTracking{ii}.position.x];
                y = [y; tempTracking{ii}.position.y];
                v = [v; tempTracking{ii}.v];
                xy_angle = [xy_angle; tempTracking{ii}.angle];
                timestamps = [timestamps; tempTracking{ii}.timestamps];
                folder{ii} = tempTracking{ii}.folder;
                samplingRate = [samplingRate; tempTracking{ii}.samplingRate];
                description{ii} = tempTracking{ii}.description;
            end
            
            % save data to ouput structure
            tracking_struct.position.x = x;
            tracking_struct.position.y = y;
            tracking_struct.v = v; 
            tracking_struct.angle = xy_angle;
            tracking_struct.folders = folder;
            tracking_struct.samplingRate = samplingRate;
            tracking_struct.timestamps = timestamps;
            tracking_struct.description = description;
            tracking_struct.vidnames = vidnames;
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
            p.addParameter('likelihood',.95,@isnumeric); % tracking quality thres [0-1]
            p.addParameter('pulses_delta_range',0.01,@isnumeric); % range for ttls
            
            p.parse(varargin{:});
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
            
            file_idx = find(contains(dlc_files,'DLC'));
            dlc_files = dlc_files(file_idx);
            vidnames = vidnames(file_idx);
            ep_idx = ep_idx(file_idx);
            frame_rate = frame_rate(file_idx);
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
                disp(['N number of frames off: ', num2str(basler_intan_diff)])
                warning('Unnoticed problem with Camera/Intan... I would go back and check both step by step');
            end
        end
        
        function [tracking_files,vid_names,ep_idx,frame_rate] =  get_tracking_info_from_session(session,basename)
            % extracts tracking information (tracking filename, video name, and epoch index
            % basename.session.behavioralTracking
            
            disp(['Checking for tracking files in ', basename, '.session.behavioralTracking'])
            
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
        
        % Utils 
        function rescale_maze_coords(basepath)
            % rescale_maze_coords converts cooridnates maze_coords.csv to
            % scaled coords using the scale params. Assumes tracking has
            % been processed process_tracking. 
            %
            %
           
             % check for maze_coords and if doesn't exist return
            if isempty(dir(fullfile(basepath,'*_maze_coords.csv')))
                warning('No maze coords found. Run get_maze_XY.m to obtain maze coords')
                return
            end
            
            % load files
            basename = basenameFromBasepath(basepath);
            session = loadSession(basepath,basename);
            load(fullfile(basepath,[basename,'.animal.behavior.mat']))
            
            % get the scale parameters for each epoch
            [x_origin,y_origin,scale_factor_x,scale_factor_y] = tracking.get_scale_params(session,behavior,basepath);

            % loop through behavioral tracking epochs, rescale xy and save
            % back to maze_coord_df. 
            for i = 1:length(session.behavioralTracking)
                % load coords from struct or csv in session folder
                if isfield(session.behavioralTracking{1, i},'maze_coords')
                    maze_coord_temp = session.behavioralTracking{1, i}.maze_coords;
                else
                    maze_coord_temp = readtable(fullfile(basepath,...
                    [extractBefore(session.behavioralTracking{1,i}.notes,'.avi'),'_maze_coords.csv']));
                end
                maze_coord_temp.x_scaled = tracking.scale_transform(maze_coord_temp.x,x_origin{i},scale_factor_x{i});
                maze_coord_temp.y_scaled = tracking.scale_transform(maze_coord_temp.y,y_origin{i},scale_factor_y{i});
                session.behavioralTracking{1, i}.maze_coords = maze_coord_temp;
                writetable(maze_coord_temp,fullfile(basepath,...
                    [extractBefore(session.behavioralTracking{1,i}.notes,'.avi'),'_maze_coords.csv']));
            end
            % save session file
            save(fullfile(basepath,[basename,'.session.mat']),'session')
        end
    end
end