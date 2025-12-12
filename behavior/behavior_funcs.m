classdef behavior_funcs
    %behavior consists of functions used to assess locomotion of rodent
    % behavior
    %
    % Adapted from OF from ephys_tools by LB 12/2022
    
    
    
    properties
        Property
    end
    
    methods(Static)
        
        %% General locomotion
        
        function length = path_length(x,y,velocity,varargin)
            % Computes total path length during running
            % inputs:
            %   x: position vector of x-coordinates of length n
            %   y: position vecotr of y-coordinates of length n
            %   velocity: vector of instantaneous velocity (length of n - 1)
            % output:
            %   length: total path length in cm for active motion
            
            p = inputParser;
            addParameter(p,'run_threshold',3,@isnumeric);
            
            parse(p,varargin{:});
            run_threshold = p.Results.run_threshold;
            
            %distance formula
            distance_vector = sqrt((diff(x)).^2 + (diff(y)).^2);
            
            % length of path when animal is running
            length = nansum(distance_vector(velocity >= run_threshold,1));%Total Path length for points greater than 3cm/s
            
        end
        
        function [stopIdx,startStop,endStop] = stop_epochs(velocity,varargin)
            % Given a vector of velocity, find when the vector falls below
            % input:
            % run_threshold: velocity that separates movement from
            %   non-movement in same unites as xy coordinates (ie cm/s or
            %   pixels/s) default: 3 cm/s
            % epoch: number of contiguous frames to search for stops.
            % Default is 60 frames assuming 60Hz sample rate (i.e. animal
            %   is stopped for 1 second).
            % output:
            %   stopIdx: logical index of time associated with stop
            %   startStop: row index for start of stops
            %   endStop: row index for end of stops
            
            p = inputParser;
            addParameter(p,'run_threshold',3,@isnumeric);
            addParameter(p,'epoch',60,@isnumeric);
            
            parse(p,varargin{:});
            run_threshold = p.Results.run_threshold;
            epoch = p.Results.epoch;
            
            stopIdx = contiguousframes(velocity <= run_threshold,epoch);
            [startStop,endStop,~] = find_groups(stopIdx);
        end
        
        function center = center_of_mass(x,y)
            %Find center of mass for set of x,y coordinates
            % uses min_encl_ellipsoid to find center of mass of coordinates
            % and results [x,y] as 'center'
            temp=[x,y];
            % remove nans
            temp(any(isnan(temp),2),:)=[];
            % min_encl_ellipsoid requires size > 1
            if ~isempty(temp) && size(temp,1)>1
                [ ~,~, center] = min_encl_ellipsoid(temp(:,1),temp(:,2));
            else
                center = [NaN,NaN];
            end
        end
        
        function stop_measures = stops(x,y,ts,velocity,fr,varargin)
            % Finds when animal stop ( < 3cm/s velocity for at least 1
            % second) and computes stop features.
            % inputs:
            %       x: position vector for x of length n
            %       y: position vector for y of length n
            %       ts: timestamps of length n
            %       fr: frame rate
            %       velocity: instantaneous velocity of length n-1
            %       epoch (optional) : length of minimum stop epoch in
            %       frames (default is 1 second (frame rate in Hz))
            % outputs:
            %       stop_measures: structure containing outcome measures.
            %           stopIdx: logical index of points where rat is
            %                    stopped.
            %           stops: cell array of stop xy coordinates
            %           timeStopped: time (seconds) of a stop
            %           tsStop: cell array of time stamps corresponding to stop
            %           NumStops: number of stops made
            %
            p = inputParser;
            addParameter(p,'run_threshold',3,@isnumeric);
            addParameter(p,'epoch',[],@isnumeric);

            
            parse(p,varargin{:});
            run_threshold = p.Results.run_threshold;
            epoch = p.Results.epoch;
            
            if isempty(epoch)
                epoch = fr;
            end

            [stop_measures.stopIdx,startStop,endStop] = behavior_funcs.stop_epochs(velocity,'epoch',epoch,'run_threshold',run_threshold);
            
            if isempty(startStop)
                disp('No stops detected, check speed threshold')
                stop_measures.NumStops = 0;
                return
            end
            
            %     This finds coordinates associated with stop
            for ii=1:length(startStop)
                motionless{ii} = [x(startStop(ii):endStop(ii)),y(startStop(ii):endStop(ii))];
                timeMotionless{ii} = size(motionless{ii},1)/fr;
                tsStop{ii} = ts(startStop(ii):endStop(ii));
                meanTsStop(ii) = mean(ts(startStop(ii):endStop(ii)));
            end
                        
            stop_measures.stops = motionless;
            stop_measures.timeStopped = timeMotionless;
            stop_measures.tsStop = tsStop;
            stop_measures.meanTsStop = meanTsStop; 
            %Create Number of Stops
            stop_measures.NumStops = size(stop_measures.stops,2);
            
            %Find center of mass for stops.
            for ii=1:size(stop_measures.stops,2)
                center(ii,:) = behavior_funcs.center_of_mass(stop_measures.stops{1,ii}(:,1),stop_measures.stops{1,ii}(:,2));
            end
            
            stop_measures.stopCenter = center;
            
        end
        
        function quad_measures = quadrant(x,y,fr,varargin)
            % Finds when animal stops ( < 3cm/s velocity for at least 1
            % second) and computes stop features.
            % inputs:
            %       x: vector of cordinates for x
            %       y: vector of cordinates for y
            %       fr: video sampling rate in Hz
            % outputs:
            %       quad_measures: structure containing outcome measures.
            %           stopIdx: logical index of points where rat is
            %                    stopped.
            %           stops: cell array of stop xy coordinates
            %           timeStopped: time (seconds) of a stop
            %           tsStop: cell array of time stamps corresponding to stop
            %           NumStops: number of stops made
            %
            % LB 2022
            p = inputParser;
            p.addParameter('n_quadrants',4,@isnumeric)
            p.addParameter('run_threshold',3,@isnumeric)
            
            p.parse(varargin{:})
            n_quadrants = p.Results.n_quadrants;
            run_threshold = p.Results.run_threshold;
            
            %squared distance between consecutive points
            sqrXDiff = (diff(x)).^2;
            sqrYDiff = (diff(y)).^2;
            
            %distance formula
            distance_vector = sqrt(sqrXDiff + sqrYDiff);
            %instanteous velocity
            velocity = (distance_vector)*fr; %instanteous velocity
            velocity = smoothdata(velocity,'sgolay',fr*.8); %smoothed with moving median window over less than 1 second
            
            clear sqrXDiff sqrYDiff
            
            [stopIdx,startStop,endStop] = behavior_funcs.stop_epochs(velocity,epoch,'run_threshold',run_threshold);
            
            % animal is stopped for in quadrant for 1 second
            stopIdx = contiguousframes(velocity <= run_threshold,fr);
            [startStop,endStop,~] = find_groups(stopIdx);
            
            %     This finds coords for stopping
            for ii=1:length(startStop)
                motionless{ii}=[back(startStop(ii):endStop(ii),1),back(startStop(ii):endStop(ii),2)];
            end
            
            stops = motionless;
            
            clear startStop endStop
            
            %Find center of mass for stops.
            for ii=1:size(stops,2)
                temp=[stops{1,ii}(:,1),stops{1,ii}(:,2)];
                temp(any(isnan(temp),2),:)=[];
                if ~isempty(temp) && size(temp,1)>1
                    [ ~,~, stopCenter] = min_encl_ellipsoid(temp(:,1),temp(:,2));
                else
                    stopCenter(1,1)=NaN;
                    stopCenter(2,1)=NaN;
                end
                stop(ii,1)=stopCenter(1,1);
                stop(ii,2)=stopCenter(2,1);
            end
            
            stopCenter=stop;
            
            clear stop
            
            %calculate verticies for quadrants
            quadrants = createZones([0,0],params.dia{j},'numQuad',n_quadrants,'fig',0); %16 pie shaped segments
            
            
            %Calculate dwell time per quadrant
            for i=1:size(quadrants,1)-1
                tempXin = [0 quadrants(i,1) quadrants(i+1,1)];
                tempYin = [0 quadrants(i,2) quadrants(i+1,2)];
                [in,~] = inpolygon(back(:,1),back(:,2),tempXin,tempYin);
                quad_measures.dwellQuad{j}(1,i) = sum(in)/fr; %
                quad_measures.pathLQuad{j}(1,i) = sum(distance_vector(in(1:end-1) & velocity>3 ,1)); %mean distance when animal is in quadrant and moving >3cm/s
                quad_measures.pathIVQuad{j}(1,i) = nanmean(velocity(in(1:end-1),1)); %mean linear velocity
                quad_measures.numstopQuad{j}(1,i) = sum(inpolygon(stopCenter{j}(:,1),stopCenter{j}(:,2),tempXin,tempYin));
                quad_measures.stopQuad{j}(1,i) = sum(in(1:end-1,:) & stopIdx{j})/fr; %stop duration
                quad_measures.angVelQuad{j}(1,i) = nanmean(abs(angVel(in(1:end-1,:))));
                clear tempXin tempYin in
            end
            
        end
        
        function [occ,map] = occ_map(x,y,binsize,fr)
            % Builds an occupancy map for maze. Assums params table.
            % inputs:
            %   x: position vector of x-coordinates of length n
            %   y: position vecotr of y-coordinates of length n
            %   binsize: how big each bin should be in cm
            %   fr: frame rate
            %
            % output:
            %   occ: smoothed occupancy map
            %   map: raw occupancy map
            
            %Creates bin edges for heatmap
            
            xedge=linspace(nanmin(x),nanmax(x),round(range(x)/binsize));
            yedge=linspace(nanmin(y),nanmax(y),round(range(y)/binsize));
            
            %Bin coordinates for heat map and apply guassian filter to smooth
            [map] = histcounts2(x,y,xedge,yedge);
            map=flipud(imrotate(map,90));
            map=map/fr;
            
            % smooths over 1.5 cm
            occ = imgaussfilt(map, 1);
            occ(map==0) = nan;
            map(map==0) = nan;
        end
        
        function [out] = thigmotaxis_square(x,y,fr, basepath, epoch)
            
            [xOuter, yOuter, xInner, yInner] = generateNestedSquares(basepath, epoch); 
             %Calculate dwell time for outter edge
            [in,~] = inpolygon(x,y,xInner,yInner);
            out = sum(~in)/fr;
            
        end
        
        function [out] = thigmotaxis_circle(x,y,fr,diameter,varargin)
            % NOT FUNCTION YET NEED TO ADJUST CREATEZONES FOR SQUARE
            % ENVIORNMENT LB 12/2022
            
            % Computes time spent near outside wall
            p = inputParser;
            addParameter(p,'center_proportion',.8,@ischar);
            
            parse(p,varargin{:});
            center_proportion = p.Results.center_proportion;
            
            %Create annuli for center and outter edge of maze
            outsideDwell = createZones([0,0],diameter,'type','annulus','fig',0,'annulusSize',center_proportion); %default annulus size is 80% for createZones
            
            %Calculate dwell time for outter edge
            [in,~] = inpolygon(x,y,outsideDwell(:,1),outsideDwell(:,2));
            out = sum(~in)/fr;
            
            
        end
        
        function sa = search_area(map,varargin)
            % computes proportion of maze area occupied by animal given
            % occupancy map.
            % inputs:
            %   map: raw occpancy map (n x n grid of binned occupancy
            %   weighted by frame rate)
            
            % output:
            %   sa: search area as a proportion (e.g. .10 = 10%
            %   of maze searched).
            p = inputParser;
            addParameter(p,'maze_type','square',@ischar);
            
            parse(p,varargin{:});
            maze_type = p.Results.maze_type;
            
            %Use meshgrid to serve as basis for logical mask.
            imageSizeX = size(map,1);
            imageSizeY = size(map,2);
            [columnsInImage,rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
            
            clear imageSizeX imageSizeY
            
            % Next create the circle in the image.
            if ismember(maze_type,'circle')
                centerX = median(1:size(map,1));
                centerY = median(1:size(map,2));
                radius = median(1:size(map,2));
                circlePixels = (rowsInImage - centerY).^2 ...
                    + (columnsInImage - centerX).^2 <= radius.^2;
                map(~circlePixels)=NaN; %indicate area outside the maze by labelling with NaN
            end
            
            clear rowsInImage columnsInImage
            
            [length, width] = size(map);
            sa = sum(sum(map>.100))/(length*width); %Calculate the proportion of bins occupied by animal.
            
        end
        
        function ED = egocentric_direction(x,y, object_center_x,object_center_y)
            
            ED = mod(atan2d(y-object_center_y,x-object_center_x),rad2deg(2*pi));

        end

        function CD = cue_direction(ED, HD)
            % takes in egocentric direction and head direction
            % (analogSignalArray objects) and computes direction to cue by
            % taking the circular distance between the two headings. 
            % returns analougSignalArray with cue direction. 
            r = angle(exp(1i*deg2rad(ED.data(:,1)))./exp(1i*deg2rad(HD.data(:,1))));
            
            % shift all by 45 
            r = angle(exp(1i*r)./exp(1i*deg2rad(45)));
            CD = rad2deg(r');
            

        end
        
        
        %% Cheeseboard 

        function process_napari_cheeseboard(basepath)
            %PROCESS_NAPARI_CHEESEBOARD Process and sync Napari-points output for cheeseboard task with SNLAB general behavior file.
            % Requires that SNLab process_tracking has been completed for
            % the basepath.
            %
            % 
            % Parameters
            % ----------
            % basepath : string
            %     Path to the base directory containing session files and
            %     CSV annotations. CSVs must be in basepath and labelled
            %     *e_rewards.csv, *e_start.csv, *_trials.csv, 
            %
            % Description
            % -----------
            % This function performs three key operations for cheeseboard behavioral data:
            % 1. Scales XY coordinates for rewards and start points drawn in Napari.
            % 2. Updates session struct with scaled coordinates.
            % 3. Syncs Napari trial annotations to behavior timestamps.
            % 
            % Note: Naming convention for filenames is the
            % video_name_*.csv, i.e.
            % 4456_learning_10072025-095012_rewards.csv
            % 
            % LB 10/2025


            
                
            % Scale XY coordinates of start, save to session and basepath
            behavior_funcs.scale_napari_coords(basepath,'csv_tag', '_rewards')
            
            % Scale XY coordinates of rewards, save to session and basepath
            behavior_funcs.scale_napari_coords(basepath,'csv_tag','_start')
            
            % Sync trials to behavior file 
            behavior_funcs.napari_trials_to_behavior(basepath)
            
        end

        function scale_napari_coords(basepath,varargin)
        %SCALE_NAPARI_COORDS Scale XY coordinates from Napari to behavioral coordinates.
        %
        % Parameters
        % ----------
        % basepath : string
        %     Base directory path containing the session and CSV files.
        %
        % Optional Parameters
        % -------------------
        % 'csv_tag' : string, default = '_rewards'
        %     Tag used to find the appropriate CSV files (e.g., '_start', '_rewards').
        %
        % Description
        % -----------
        % This function reads CSV annotations made in Napari, scales the XY pixel
        % coordinates based on the reference distance from the session struct, and
        % stores the scaled data back into the session struct and as new CSV files.
        % LB 10/2025

            % input parser
            p = inputParser;
            p.addParameter('csv_tag','_rewards',@ischar);
            
            
            p.parse(varargin{:}) 
            csv_tag = p.Results.csv_tag; 
            
            
            basename = basenameFromBasepath(basepath);
            session = loadSession(basepath,basename);
            csv_files = my_dir(fullfile(basepath,['*',csv_tag,'.csv']));

            if isempty(csv_files)
                error(['no files found in ', basepath, ': with ',csv_tag,' included in the file extension'])
            end

            
            % loop through video
            for file = 1:length(session.behavioralTracking) %loop through folders containing subject videos
               
                vidname = extractBefore(session.behavioralTracking{1,file}.notes,'.'); % default video name and file extension saved in notes for SnLab process tracking 
                reward_table = readtable(fullfile(basepath,csv_files.name{contains(csv_files.name,vidname)}));
            
                % load pixel distance and pixel_reference 
                pixel_distance = session.behavioralTracking{1, file}.pixel_distance;
                pixel_dist_reference = session.behavioralTracking{1, file}.pixel_dist_reference;
            
                reward_table.x_scaled = reward_table.axis_2 * (pixel_dist_reference/pixel_distance);
                reward_table.y_scaled = reward_table.axis_1 * (pixel_dist_reference/pixel_distance);
                
                % save to session 
                session.behavioralTracking{1,file}.(['cheeseboard',csv_tag]) = reward_table;
                
                % save data to csv 
                save_file = fullfile(basepath,[vidname,csv_tag,'.csv']);
                writetable(reward_table,save_file);
            end
            
            % save session back to basepath 
            save(fullfile(basepath,[basename,'.session.mat']),'session');
        
        end
        
        function napari_trials_to_behavior(basepath)
            %NAPARI_TRIALS_TO_BEHAVIOR Convert Napari trial CSVs (points layer) into CellExplorer animal behavior file format.
            %
            % Parameters
            % ----------
            % basepath : string
            %     Path to the session base directory, containing csvs,
            %     animal.behavior.mat file, session.mat file
            %
            % Description
            % -----------
            % This function loads trial annotations (e.g., start and stop frames) from Napari-drawn
            % CSV (*_trials.csv) and syncs them with xy coordinates saved in the animal behavior`.mat` file.
            % It updates `behavior.trials` and `behavior.trialsID` with synchronized timestamps and epoch_ID_trialN, respectively.
            %
            % LB 10/05/2025


            % crate basename from basepath
            basename = basenameFromBasepath(basepath);
            
            % load session and animal behavior file
            load(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior');
            
            % load session file for behavioralTracking that contains filenames
            session = loadSession(basepath,basename);
            epoch_df = load_epoch('basepath', basepath);
            
            % find all csv files that contain *_rewards,*_trials,*_start
            csv_files = my_dir(fullfile(basepath,'*_*.csv'));
            csv_files = csv_files(contains(csv_files.name,{'_rewards','_start','_trials'}),:);
            
            
            trials_mat = [];
            trials_id = {};
            % loop through videos indicated in session.behavioralTracking
            for ii = 1:length(session.behavioralTracking)
            
                % load epoch to get start/end for timestamps
                epoch = session.behavioralTracking{1,ii}.epoch;
                name = session.epochs{1,epoch}.name;
                % grap video so we can find row index for a given video
                vidname = session.behavioralTracking{1,ii}.notes;
                vidname = extractBefore(vidname,'.');
                
                trials_idx = contains(csv_files.name,vidname) & contains(csv_files.name,"trials");
            
                % load epoch trials
                trials_df = readtable(fullfile(basepath,csv_files.name{trials_idx}));
                
                [~, idx] = sort(trials_df.axis_0);
                trials_df = trials_df(idx, :);
                
                idx = InIntervals(behavior.timestamps, [epoch_df.startTime(epoch), epoch_df.stopTime(epoch)]);
                current_ts = behavior.timestamps(idx);
                starts = (trials_df.axis_0(1:2:end) / behavior.sr) + current_ts(1);
                stops = (trials_df.axis_0(2:2:end) / behavior.sr) + current_ts(1);
               
                trials_mat = [trials_mat; starts stops];
                temp_id = {};
            
                for t = 1:length(starts)
                    temp_id{t} = [name,'_trial',num2str(t)];
                end
            
                trials_id = [trials_id; temp_id'];
            
            end
    
            % must be equal dims to add
            assert(size(trials_mat,1) == length(trials_id))
            
            % if trials exist, which they should do if behavior file dose then
            % concatenate
            if ~isempty(behavior.trials)
                % load and add new
                merged_trial_mat = [behavior.trials; trials_mat];
                merged_id = [behavior.trialsID; trials_id];
                
                % overwrite with added trials
                behavior.trials = merged_trial_mat;
                behavior.trialsID = merged_id;
            else
                % Add to behavior file
                behavior.trials = trials_mat;
                behavior.trialsID = trials_id;
            end
            % save behavior file
            save(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')
        
        end

        %% Home base functions
        
        function metrics = home_base_metics(out_home,x,y,velocity,x_home,y_home,stopIdx,tsStop)
            
            % This finds stops that occur in the home base boundary
            [startStop,~,~] = find_groups(stopIdx);
            
            [~,~,metrics.entries] = find_groups(out_home);
            
            % total time in home base
            metrics.hbOcc = nansum(out_home)/fr;
            
            % average velocity in home base
            metrics.hbVel = nanmean(velocity(out_home(1:end-1,1),1)); %remove last tempIn idx to accomodate velocity length
            
            
            % time moving slow in home base in seconds
            metrics.slowInHB = nansum(out_home(1:end-1) & stopIdx)/fr; % time being slow in homebase
            
            % proportion of time being slow in home base
            metrics.HBclass= slowInHB/hbOcc; % proportion of time being slow in hb
            
            % number times an animal stoped in the home base
            metrics.HBstops = nansum(inpolygon(x(startStop,1),...
                y(startStop,2),x_home(1:end-2)',y_home(1:end-2)')); % Find number of times animals initiated a start in the home base
            
            % index for the time it takes to reach the home base
            time2HB_idx = inpolygon(x(startStop,1),...
                x(startStop,2),x_home(1:end-2)',y_home(1:end-2)');
            
            % time to first time in home base
            firstIdx = find(time2HB_idx);
            
            % save time to home base
            time2HB_time = tsStop{1,firstIdx(1)};
            metrics.time2HB = time2HB_time(1,1);
            
            
        end
        
        function [HBBound,HBcoords,HBcenter,out_home,x_home,y_home] = rescale_home_base(home_base_x,home_base_y,upHBmap,upFactor,diameter,fr)
            
            %rescale coordinates back to pool size
            x_home = rescale([home_base_x,upFactor+1,size(upHBmap,2)-upFactor],-(diameter/2),(diameter/2));
            y_home = rescale([home_base_y,upFactor+1,size(upHBmap,2)-upFactor],-(diameter/2),(diameter/2));
            
            % logical index for all coordinates inside home base
            in_home = inpolygon(x,y,x_home(1:end-2)',y_home(1:end-2)');
            
            % coordinates for frames inside home base for at least 2seconds
            out_home = contiguousframes(in_home,fr*2); %has to be inside of hb for at least 2 sec to count as entryd
            
            % boundary of home base
            HBBound = [x_home(1:end-2)',y_home(1:end-2)'];
            
            % time in home base
            HBcoords = [x(out_home,1),y(out_home,2)];
            
            [ ~,~, tempC] = min_encl_ellipsoid(HBBound(:,1),HBBound(:,2));
            HBcenter= [tempC(1,1),tempC(2,1)];
        end
        
        function [HB_stop_dist,HB_close_stop,HB_stop_dist_vector] = stops_to_homebase(HBcenter,stops)
            
            % Average Proximity of stops from hb center
            disttemp = zeros(size(stops,2),1);
            for i = 1:size(stops,2)
                temp = stops{1,i};
                dist = sqrt((HBcenter(1,1)-temp(:,1)).^2+(HBcenter(1,2)-temp(:,2)).^2);
                disttemp(i,1) = nanmean(dist);
                
            end
            
            HB_stop_dist = nanmean(disttemp); %average stop distance
            HB_close_stop = sum(disttemp<25); %number of stops within 25cm from hb center
            HB_stop_dist_vector = disttemp;
            
        end
        
        function [HB_avg_dist,HB_max_dist,HB_min_dist] = distance_between_homebase(HBcenter)
            % input:
            %   - HBcenter: vector of HB centers
            % output:
            %   -HB_avg_dist: average distance between home bases
            %   -HB_max_dist: maximum distance between home bases
            %   -HB_min_dist: minimum distance between home bases
            %
            % Calculate distance measures between high occupancy coordinates centers
            temp = [];
            if size(HBcenter,2) > 1
                for hb = 1:size(pHBcenter,2)
                    temp = [temp; HBcenter{1,hb}];
                end
                HB_avg_dist = nanmean(pdist(temp));
                HB_max_dist = max(pdist(temp));
                HB_min_dist = min(pdist(temp));
            else
                HB_avg_dist = NaN;
                HB_max_dist = NaN;
                HB_min_dist = NaN;
            end
        end

        function HBdist2Cue =  homebase_dist_from_cue(cueCM,HBcenter)
            % calculating proximity of high occupancy coordinates center from the
            % cue boundary
            
            k = convhull(cueCM(:,1),cueCM(:,2));
            cueBoundary = [cueCM(k,1),cueCM(k,2)];
            
            for r = 1:size(HBcenter,2)
                distances = sqrt(sum(bsxfun(@minus, cueBoundary, [HBcenter{1,r}(:,1),HBcenter{1,r}(:,2)]).^2,2));
                HBdist2Cue{1,r} = unique(distances(distances == min(distances))); %Find minimum distance from cue boundary to hb center
            end
            
        end
        
        
        %% Object exploration function
        function [results,results_binned,explore_vectors] = score_object_exploration(basepath,varargin)
            % within basepath, scores epochs with 'object' indicated in
            % session.epoch.name. Determines the time spent exploring
            % objects by positions listed in
            % session.behavioralTracking.maze_coords.
            %
            % Animal scoring adapted from Leger et al 2013, with exception
            % that climbing is included  as exploration.
            % https://doi.org/10.1038/nprot.2013.155
            
            % Output:
            %   structure containing,
            %   results: first object explored (string),
            %   time to object exploration (s), discrimination ratio for
            %   first 180 and 300 seconds after first object is explored.
            
            %   vectors: object A & B min distance to object boundary.
            
            % Assumes general behavior file and *maze_coords.csv is in basepath.
            p = inputParser;
            p.addParameter('distance_threshold',6,@isnumeric) %4cm as tracking is of back of electrode cage for ephys experiments
            p.addParameter('duration_in_seconds',90,@isnumeric)
            p.addParameter('heading_thresh',360,@isnumeric)
            p.addParameter('bin_width',90,@isnumeric)% in seconds
            p.addParameter('fig',true,@islogical)
            
            p.parse(varargin{:});
            distance_threshold = p.Results.distance_threshold;
            duration_in_seconds = p.Results.duration_in_seconds;
            heading_thresh = p.Results.heading_thresh;
            bin_width = p.Results.bin_width;
            fig = p.Results.fig;

           
          
            [~, basename] = fileparts(basepath);
            
            % load animal behavior and session files
            session = loadSession(basepath);
            try
                load(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')
            catch
                disp(['No behavior file found for basepath ',basepath,'. Exiting function'])
                results = NaN;
                explore_vectors = NaN;
                return
            end
            
            % load epochs
            behave_ep = behavior_funcs.load_epochs(basepath);
            trial_ep = behavior_funcs.load_trials(basepath);
            
            % get xy as analog signal array for easy epoching
            positions = analogSignalArray(...
                'data',[behavior.position.x;behavior.position.y],...
                'timestamps',behavior.timestamps,...
                'sampling_rate',behavior.sr);
            
            % load animal head direction
            HD = behavior_funcs.load_HD(basepath);
  
            % determine which object moved between the two epochs 
            moved_object_id =  behavior_funcs.find_moved_object(basepath);
            
            % initialize output 
            results = table;
            results.basepath = repmat(basepath,behave_ep.n_intervals,1);
            results_binned = cell2table(cell(1, 7), 'VariableNames', {'bin_n','bin_center_second','object_id','object_explore','epoch_name','basepath','object_identity'});
            % for each behavior epoch, get exploration of objects within
            % epoch
            for ep = 1:behave_ep.n_intervals
                
                epoch_name{ep} = session.epochs{1,session.behavioralTracking{1, ep}.epoch}.name;
                
                % get intervalArray of current epoch
                cur_ep = behave_ep(ep) & trial_ep;
                
                % animal position in epoch
                cur_pos = positions(cur_ep);
                cur_HD = HD(cur_ep);
                
                            % get maze coords for object coordinates
                [object_center,object_edge,object_id] = behavior_funcs.load_object_coords_for_epoch(session,behavior,ep);
                
                % compute cue angle for each object 
                for obj = 1:size(object_center,1)
                    center_temp = object_center(obj,:);
                    cur_ED = behavior_funcs.egocentric_direction(cur_pos.data(:,1),cur_pos.data(:,2), center_temp(1),center_temp(2));
                                        % get xy as analog signal array for easy epoching
                    cur_ED = analogSignalArray(...
                        'data',cur_ED,...
                        'timestamps',cur_pos.timestamps,...
                        'sampling_rate',cur_pos.sampling_rate);
                    cue_angle(obj,:) = behavior_funcs.cue_direction(cur_ED, cur_HD);
                end

                % cue angle 
                % get object boundary
                [bound_x, bound_y] = behavior_funcs.object_boudary(object_center,object_edge);
                
                % get distances from objects
                object_distances =  behavior_funcs.dist_from_object(cur_pos.data(:,1)...
                    ,cur_pos.data(:,2)...
                    ,bound_x,bound_y);
                
                obj_idx = contains(object_id,'A');

                % get exploration time for each object 
                [explore_time,explore_vec] =  behavior_funcs.object_explore(cur_pos,...
                        cue_angle,...
                        heading_thresh,...
                        object_distances,...
                        distance_threshold);
                
                
                % get the time for first object exploration 
                 [start_explore, idx] = min([min(explore_vec{obj_idx}), min(explore_vec{~obj_idx})]);
                 
                % Get the object they first visited, and time it took to
                % visit it
                if ~isempty(idx) 
                    first_object_explored{ep} = object_id{idx};
                    time_2_object_explore(ep) = start_explore - cur_pos.timestamps(1); % first object explore ts minus start timestamps

                else % if they didn't explore any objects, initalize all the output variables with nan
                    
                    first_object_explored{ep} = 'None';
                    time_2_object_explore(ep) = NaN;
                    object_explore_restrict(1,:) = NaN;
                    object_explore_restrict(2,:) = NaN;
                    
                    obj_A_explore(ep) = NaN;
                    obj_B_explore(ep) = NaN;
                    obj_A_explore_restrict(ep) = NaN;
                    obj_B_explore_restrict(ep) = oNaN;

                   % store vectors for each epoch 
                   explore_vectors{ep}.task_id = session.epochs{1,session.behavioralTracking{1, ep}.epoch}.name;
                   explore_vectors{ep}.object_id = object_id;
                   explore_vectors{ep}.bin_explore = NaN;
                   explore_vectors{ep}.bin_explore_table = NaN;
                   explore_vectors{ep}.explore_vec = NaN;
                   explore_vectors{ep}.cue_angle = NaN;
                   explore_vectors{ep}.object_distances = NaN;
                   explore_vectors{ep}.start_explore = NaN;
               
                    bin_explore_table = cell2table(cell(1, 7), 'VariableNames', {'bin_n','bin_center_second','object_id','object_explore','epoch_name','basepath','object_identity'}); 
                    bin_explore_table.epoch_name =  repmat(epoch_name(ep),size(bin_explore_table,1),1);
                    bin_explore_table.basepath =  repmat(basepath,size(bin_explore_table,1),1);
                    bin_explore_table.object_identity =  repmat('fixed',size(bin_explore_table,1),1);
                    bin_explore_table(ismember(bin_explore_table.object_id, moved_object_id),:).object_identity = repmat('moved',sum(ismember(bin_explore_table.object_id, moved_object_id)),1);
                    
                    clear cue_angle object_distances bin_explore explore_vec
                    
                    continue

                end
                  
                
                for obj = 1:length(object_id)
                    object_explore_restrict(obj,:) = behavior_funcs.limit_explore_to_segment(start_explore,...
                        duration_in_seconds,explore_vec{obj});
                end
                 
                
                % limit position for binned data to include only time after
                % animal started exploring objects 
                explore_epoch = IntervalArray([start_explore,cur_pos.timestamps(end)]);
                
                % put object distances and cue angle in array 
                 cue_angle = analogSignalArray(...
                'data',cue_angle,...
                'timestamps',cur_pos.timestamps,...
                'sampling_rate',cur_pos.sampling_rate);
            
                 object_distances = analogSignalArray(...
                'data',object_distances,...
                'timestamps',cur_pos.timestamps,...
                'sampling_rate',cur_pos.sampling_rate);
            
% %                 exploration over time
%                 [bin_explore, bin_explore_table] = behavior_funcs.object_explore_over_time(cur_pos(explore_epoch),...
%                     bin_width,... % time in seconds
%                     cue_angle(explore_epoch),... % instantaneous egocentric angle from objects
%                     heading_thresh,...% threashold for heading 
%                     object_distances(explore_epoch),...% instantaneous distance from objects
%                     distance_threshold,... % threshold for close enought 
%                     object_id); % ID for each object 
%                 
%                bin_explore_table.epoch_name =  repmat(epoch_name(ep),size(bin_explore_table,1),1);
%                bin_explore_table.basepath =  repmat(basepath,size(bin_explore_table,1),1);
%                bin_explore_table.object_identity =  repmat('fixed',size(bin_explore_table,1),1);
%                bin_explore_table(ismember(bin_explore_table.object_id, moved_object_id),:).object_identity = repmat('moved',sum(ismember(bin_explore_table.object_id, moved_object_id)),1);


%                results_binned = [results_binned; bin_explore_table];

                % Discrimination ratio for entire session, and restricted
                % session. Note - restriction is only applied to object
                % test. 
                if ismember(moved_object_id,'A') && (contains(epoch_name{ep},'object_learning'))
                    DR_overall(ep) = behavior_funcs.discrimination_ratio(explore_time(~obj_idx),explore_time(obj_idx));
                    DR_restrict(ep) = DR_overall(ep);
                    
                elseif ismember(moved_object_id,'A') && (contains(epoch_name{ep},'object_test'))
                    DR_overall(ep) = behavior_funcs.discrimination_ratio(explore_time(~obj_idx),explore_time(obj_idx));
                    DR_restrict(ep) = behavior_funcs.discrimination_ratio(object_explore_restrict(~obj_idx),object_explore_restrict(obj_idx));
                    
                elseif ismember(moved_object_id,'B') && (contains(epoch_name{ep},'object_learning'))
                    DR_overall(ep) = behavior_funcs.discrimination_ratio(explore_time(obj_idx),explore_time(~obj_idx));
                    DR_restrict(ep) = DR_overall(ep);
                    
                 elseif ismember(moved_object_id,'B') && (contains(epoch_name{ep},'object_test'))
                    DR_overall(ep) = behavior_funcs.discrimination_ratio(explore_time(obj_idx),explore_time(~obj_idx));
                    DR_restrict(ep) = behavior_funcs.discrimination_ratio(object_explore_restrict(obj_idx),object_explore_restrict(~obj_idx));
                end
                
               % save for results
                obj_A_explore(ep) = explore_time(obj_idx);
                obj_B_explore(ep) = explore_time(~obj_idx);
                obj_A_explore_restrict(ep) = object_explore_restrict(obj_idx);
                obj_B_explore_restrict(ep) = object_explore_restrict(~obj_idx);

               % store vectors for each epoch 
               explore_vectors{ep}.task_id = session.epochs{1,session.behavioralTracking{1, ep}.epoch}.name;
               explore_vectors{ep}.object_id = object_id;
%                explore_vectors{ep}.bin_explore = bin_explore;
%                explore_vectors{ep}.bin_explore_table = bin_explore_table;
               explore_vectors{ep}.explore_vec = explore_vec;
               explore_vectors{ep}.cue_angle = cue_angle.data;
               explore_vectors{ep}.object_distances = object_distances.data;
               explore_vectors{ep}.start_explore = start_explore;
               
               
               
               clear cue_angle object_distances bin_explore explore_vec
            end

            
            % store results
            results.epoch = epoch_name';
            results.moved_object_id =repmat(moved_object_id,behave_ep.n_intervals,1);
            results.first_object_explored = first_object_explored';
            results.time_2_object_explore = time_2_object_explore';
            results.DR_overall = DR_overall';
            results.DR_restrict = DR_restrict';
            results.object_A_explore = obj_A_explore';
            results.object_B_explore = obj_B_explore';
            results.object_A_explore_restrict = obj_A_explore_restrict';
            results.object_B_explore_restrict = obj_B_explore_restrict';
            
            if fig
                plot_object_explore_results(basepath,explore_vectors,results)
            end
            
        end
        
        function [explore_time_obj,explore_vec_obj] =  object_explore(cur_pos,cue_angle,heading_thresh,object_distances,distance_threshold)
            
            
            for i = 1:size(object_distances,1)
                % gather point statitics
                % time near each object (time when animal is within 5cm of object)
                [explore_time, explore_vec] = behavior_funcs.time_near_object(cur_pos,...
                    object_distances(i,:),...
                    distance_threshold,...
                    cue_angle(i,:),...
                    heading_thresh);
                
                explore_time_obj(i,:) = explore_time;
                explore_vec_obj{i,1} = explore_vec';
                
            end
        end
        
        function [bin_explore, out_table] = object_explore_over_time(ts, bin_width,cue_angle,heading_thresh,object_distances,distance_threshold,object_id)
        % FUNCTION: object_explore_over_time
        %
        % DESCRIPTION:
        % This function computes and bins object exploration within time bins, based on
        % spatial proximity and heading angle constraints. It generates a time-binned
        % summary of exploration data for multiple objects, outputting the binned
        % exploration times and a formatted table summarizing the results.
        %
        % INPUTS:
        % - ts: timestamps in seconds, vector or analouge signal array (`pos.timestamps`).
        % - bin_width: A scalar specifying the width of each time bin (in seconds).
        % - cue_angle: A matrix of cue angles (in degrees) for each object over time.
        % - heading_thresh: A scalar threshold for acceptable heading angles (in degrees).
        % - object_distances: A matrix of distances between the subject and each object over time.
        % - distance_threshold: A scalar specifying the maximum distance (in the same units as `object_distances`) 
        %                       to consider exploration of an object.
        % - object_id: A cell array of strings specifying the IDs of the objects being analyzed.
        %
        % OUTPUTS:
        % - bin_explore: A matrix where each row corresponds to an object and each column 
        %                contains the total exploration time (in seconds) within a time bin.
        % - out_table: A MATLAB table summarizing the results, containing the following columns:
        %   * `bin_n`: The bin number for each time bin.
        %   * `bin_center_second`: The center of each time bin (in seconds).
        %   * `object_id`: The ID of the object corresponding to each bin.
        %   * `object_explore`: The exploration time (in seconds) for each object in each bin.
        %
        % USAGE:
        % [bin_explore, out_table] = object_explore_over_time(pos, bin_width, cue_angle, heading_thresh, object_distances, distance_threshold, object_id);
        %
        % NOTES:
        % - The function assumes `object_distances` and `cue_angle` have the same number of rows, 
        %   corresponding to the number of objects being analyzed.
        % - Exploration times are computed using the `behavior_funcs.explore_time` function, 
        %   which should return the total exploration time for timestamps within each bin.
        %
        % EXAMPLE:
        % pos.timestamps = 0:0.1:100; % Example timestamps
        % bin_width = 10; % 10-second bins
        % cue_angle = rand(2, length(pos.timestamps)) * 180 - 90; % Random cue angles between -90 and 90 degrees
        % heading_thresh = 30; % Heading angle threshold of 30 degrees
        % object_distances = rand(2, length(pos.timestamps)) * 10; % Random distances between 0 and 10 units
        % distance_threshold = 5; % Distance threshold of 5 units
        % object_id = {'ObjectA', 'ObjectB'}; % Object IDs
        %
        % [bin_explore, out_table] = object_explore_over_time(pos, bin_width, cue_angle, heading_thresh, object_distances, distance_threshold, object_id);

        
            if (ismember(class(ts),'analogSignalArray'))

                ts = ts.timestamps;

            end 

            if (ismember(class(cue_angle),'analogSignalArray'))

                cue_angle = cue_angle.data';

            end 

            if (ismember(class(object_distances),'analogSignalArray'))

                object_distances = object_distances.data';

            end 
        
     
           
           for obj = 1:size(object_distances,1) 
               
            % create bins with
            edges = ts(1):bin_width:ts(end); % edges for binning data
            bin_centers = ((edges(1:end-1) + diff(edges)/2) - min(edges)); % bin centers in seconds 
            angle_idx = cue_angle(obj,:) >= -heading_thresh & cue_angle(obj,:) <= heading_thresh;

            explore_ts = ts((object_distances(obj,:) <= distance_threshold) & angle_idx);
             
            for i = 1:length(edges)-1
                idx = explore_ts >= edges(i) & explore_ts <= edges(i+1);
                bin_explore(obj,i) = behavior_funcs.explore_time(explore_ts(idx));
            end
           end
           
            out_table = table; 
            out_table.bin_n = [1:size(bin_explore,2), 1:size(bin_explore,2)]';
            out_table.bin_center_second = [bin_centers'; bin_centers'];
            out_table.object_id = [repmat(object_id{1},size(bin_explore,2),1); repmat(object_id{2},size(bin_explore,2),1)];
            out_table.object_explore = [bin_explore(1,:)'; bin_explore(2,:)']; 

        end
        
        function moved_object_id =  find_moved_object(basepath)
            % determines which object was moved based on distance between learning
            % and test center coordinates.
            %
            % input
            
            % load behavior epochs
            basename =  basenameFromBasepath(basepath);
            behave_ep = behavior_funcs.load_epochs(basepath);
            session = loadSession(basepath);
            load(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')
            
            
            for ep = 1:behave_ep.n_intervals
                % get maze coords for object coordinates
                [object_center{ep},~,object_id{ep}] = behavior_funcs.load_object_coords_for_epoch(session,behavior,ep);
            end
            
            % euclidean distance between center coordinates
            epoch_dist = sqrt((object_center{1}(:,1) - object_center{2}(:,1)).^2 +...
                (object_center{1}(:,2) - object_center{2}(:,2)).^2);
            
            % find max
            [~,idx] = max(epoch_dist);
            
            % identify object
            moved_object_id = object_id{1}{idx};
            
            
            
        end
        
        function [object_explore,explore_ts] = time_near_object(pos,distance_vector,dist_thresh,egocentric_heading,heading_thresh)
            % calcultes the time near an object
            % input:
            %  - pos: analogSignalArray of position (length n)
            %  - distance_vector: instantaneous distance from object (length n)
            %  - dist_thresh: threshold of distance in same units as position.
            %  - egocentric_heading: relative heading of animal from object center
            %  (degrees)
            %  - heading_thresh: value in degrees to limit (ie 45 will give +/-45
            %  from object).
            % output:
            %  object_explore: total time in units of timestamps that reflects
            
            
            angle_idx = egocentric_heading >= -heading_thresh & egocentric_heading <= heading_thresh;
            dist_idx = distance_vector <= dist_thresh;
            explore_ts = pos.timestamps(dist_idx & angle_idx);
            object_explore = behavior_funcs.explore_time(explore_ts);
            
            
        end
        
        function [bound_x, bound_y] = object_boudary(object_center,object_edge)
            % Given the center and edge of an object, compute and return xy
            % coordinates representing the radius around an object
            % input:
            % object_center: [x,y] from maze_coords.csv
            % object_edge: [x, y]
            % center of object
            
            
            % given camera position to objects, the radius may be
            % different. We'll compile the radius for both to calculate an
            % average that will be used for the boundary. 
            for i = 1:size(object_center,1)
                x0 = object_center(i,1);%%origin
                y0 = object_center(i,2);%%origin
                
                % radius (difference from center to edge
                r(i) = abs(sqrt((object_center(i,1)-object_edge(i,1))^2+(object_center(i,2)-object_edge(i,2)).^2));
                
            end
            
            % get mean of radius 
            r = mean(r);
            
            % get boundary position for each object
            for i = 1:size(object_center,1)
                x0 = object_center(i,1);%%origin
                y0 = object_center(i,2);%%origin
               
                % create circle of 25 points around object
                bound_x(i,:) = cos(linspace(-pi,pi,25)) * r + x0;
                bound_y(i,:) = sin(linspace(-pi,pi,25)) * r + y0;
                
            end
   
            
        end
        
        function metrics = object_explore_vectors(object_center,object_edge,behavior,varargin)
            % calculate object interaction for all xy coordinates and
            % return as structure.  LB 2022
            
            p = inputParser;
            p.addParameter('near_obj_dist',3,@isnumeric) % perimeter around object in units as coordinates (default cm)
            p.addParameter('run_thres',4,@isnumeric) % (cm/s)
            p.addParameter('near_object_thresh',0.5,@isnumeric) % minimum time by object to count (seconds)
            
            p.parse(varargin{:})
            near_obj_dist = p.Results.near_obj_dist;
            run_thres = p.Results.run_thres;
            near_object_thresh = p.Results.near_object_thresh;
            % get behavior epoch
            
            % grab xy coordinates
            x = behavior.position.x;
            y = behavior.position.y;
            
            % compute velocity
            [velocity, ~,~] = linear_motion(x,y,behavior.sr,1);
            
            % get object boundary
            [bound_x, bound_y] = behavior_funcs.object_boudary(object_center,object_edge);
            
            % distance of animal relative to object boundaries
            object_distances =  behavior_funcs.dist_from_object(x,y,[bound_x, bound_y]);
            
            % find distances when animal was close to object and moving
            % slow
            explore_obj_idx = (object_distances <= near_obj_dist) &  (velocity < run_thres);
            
            % logical index for all coordinates inside object boundary
            in_home = inpolygon(x,y,bound_x',bound_y');
            
            % time moving slow for at least 1 second
            stopIdx = contiguousframes(velocity < run_thres,behavior.sr);
            
            % coordinates for frames inside the object area for at least near_object_thresh seconds
            out_home = contiguousframes(in_home,behavior.sr*near_object_thresh); %has to be inside of hb for at least 1 sec to count as entryd
            
            % This finds stops that occur in the home base boundary
            [startStop,~,~] = find_groups(stopIdx);
            
            % get n entries into object area
            [~,~,metrics.entries] = find_groups(out_home);
            
            % total time near object
            metrics.obj_Occ = nansum(out_home)/behavior.sr;
            
            % average velocity near object
            metrics.obj_Vel = nanmean(velocity(out_home(1:end-1,1),1)); %remove last tempIn idx to accomodate velocity length
            
            % time moving slow near object in seconds
            metrics.obj_slow = nansum(out_home(1:end-1) & stopIdx)/behavior.sr; % time being slow in homebase
            
            % proportion of time being slow near object
            metrics.obj_class= slowInHB/obj_Occ; % proportion of time being slow in hb
            
            % number times an animal stoped in the object area
            metrics.HBstops = nansum(inpolygon(x(startStop,1),...
                y(startStop,2),bound_x(1:end-2)',bound_y(1:end-2)')); % Find number of times animals initiated a start in the home base
            
        end
        
        function dist_vec = stop_dist_from_object(object_center,stops)
            % mesaure time moving slow/stopped around object
            % input:
            %   object_center: xy coordinates of object
            %   stops: cell array containing xy coordinates of each stop.
            % output:
            %   dist_vec: vector of length n stops containing mean distance
            %   in xy units.
            
            % Average Proximity of stops from hb center
            dist_vec = zeros(size(stops,2),1);
            for i = 1:size(stops,2)
                temp = stops{1,i};
                dist = sqrt((object_center(1,1)-temp(:,1)).^2+(object_center(1,2)-temp(:,2)).^2);
                dist_vec(i,1) = nanmean(dist); % mean distance of each stop to object
                
            end
            
        end
        
        function object_distances =  dist_from_object(x,y,bound_x,bound_y)
            % calculating distance of xy to object boundary.
            
            for i = 1:size(bound_x,1)
                for ii = 1:length(bound_x(i,:))
                    distances(:,ii) = sqrt(sum(bsxfun(@minus, [bound_x(i,ii)',bound_y(i,ii)'], [x,y]).^2,2));
                end
                object_distances(i,:) = min(distances,[],2)'; %Find minimum distance from cue boundary to animal location
            end
            
        end
        
        function [object_center,object_edge,object_id] = load_object_coords_for_epoch(session,behavior,idx)
            % loads the object center and edge within a given behavior
            % epoch.
            % input:
            % session: CellExplorer session file
            % behavior: CellExplorer behavior file
            % idx: index of cell position for
            % session.behavioralTracking
            
            % output:
            % object_center: [x, y] for each object in maze_coords
            % object_edge: [x, y] for edge of obejct in maze_coords
            % object id: cell array containing string indicating object
            % (either 'A' or 'B')
            
            % get maze coords for object coordinates
            maze_coords = table;
            for i = 1:length(session.behavioralTracking)
                ep_coords = session.behavioralTracking{1, i}.maze_coords;
                ep_coords.behave_epoch = repmat(session.behavioralTracking{1, i}.epoch,size(ep_coords,1),1);
                
                maze_coords = [maze_coords; ep_coords];
            end
            maze_coords = session.behavioralTracking{1, idx}.maze_coords;
            obj_idx = contains(maze_coords.object,{'A','B'});
            obj_center = contains(maze_coords.position,{'center'});
            object_id = unique(maze_coords.object(obj_idx));
            
            % fetch xy for center of object
            if ~contains(behavior.position.units,'pixels')
                % load scaled position as scaled coordinates are not in
                % pixels
                object_center = [maze_coords.x_scaled(obj_idx & obj_center), ...
                    maze_coords.y_scaled(obj_idx & obj_center)];
                object_edge = [maze_coords.x_scaled(obj_idx & ~obj_center), ...
                    maze_coords.y_scaled(obj_idx & ~obj_center)];
            else
                warning('WARNING Behavior units are in pixels')
                %load raw coords (in pixels)
                object_center = [maze_coords.x(obj_idx & obj_center), ...
                    maze_coords.y(obj_idx & obj_center)];
                object_edge = [maze_coords.x(obj_idx & ~obj_center), ...
                    maze_coords.y(obj_idx & ~obj_center)];
            end
            
        end
        %         %% locomotion over time
        %         function locomotion_over_epochs(x,y,[startIdx,endIdx])
        %             % examine path length, search area, velocity, and stops as a
        %             % function of epoch (could be trials, behavior epochs, or time
        %             % bins
        %             disp('NOT FUNCTIONTIONAL LB 12/2022')
        %         end
        %
        
        function DR = discrimination_ratio(object_A_explore,object_B_explore)
            % returns the discrimination ratio that represents the relative
            % proportion of time spent exploring object B relative to all
            % object exploration.
            % input:
            %   object_A_explore: total time exploring object A in seconds
            %   object_B_explore: total time exploring object B in seconds
            
            DR = (object_B_explore - object_A_explore) / (object_B_explore + object_A_explore);
            
            
        end
        
        function object_explore = limit_explore_to_segment(start_time,duration_in_seconds,object_explore_vec)
            % limits explore vector to a specific segment, for example for the
            % first five minutes after object exploration.
            object_explore_vec = object_explore_vec(object_explore_vec >= start_time & object_explore_vec <= (start_time+duration_in_seconds));
            
            object_explore = behavior_funcs.explore_time(object_explore_vec);
            
        end
        
        function object_explore = explore_time(explore_ts)        

            % for a set of timestamps corresponding to object exploration
            % 'explore_ts', explore time computes the overall time and removes
            % jumps by keeping only frames with frame rates that match the mode.
            
            delta_explore = diff(explore_ts);
            delta_explore(delta_explore > mode(delta_explore)) = nan;
            object_explore = nansum(delta_explore);
        end
        
        
        %%  behavior utils

        function trial_ep = load_trials(basepath, truncate_time)
            basename = basenameFromBasepath(basepath);
            load(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior');
            
            if isemtpy(truncate_time)
            
                trial_ep = IntervalArray(behavior.trials);
            else
                behavior.trials
            end
        end
        
        function behave_ep = load_epochs(basepath)
            % load session
            session = loadSession(basepath);
            
            % loop through behavioral epochs
            behave_ep = [];
            for ep = 1:length(session.behavioralTracking)
                behave_ep(ep,:) = [session.epochs{1,session.behavioralTracking{1,ep}.epoch}.startTime, ...
                    session.epochs{1,session.behavioralTracking{1,ep}.epoch}.stopTime];
            end
            behave_ep = IntervalArray(behave_ep);
            
        end
        
        function heading_direction = load_HD(basepath)
            % load behavior
            basename = basenameFromBasepath(basepath);
            load(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')
            % get xy as analog signal array for easy epoching
            heading_direction = analogSignalArray(...
                'data',rad2deg(behavior.angle),...
                'timestamps',behavior.timestamps,...
                'sampling_rate',behavior.sr);
            
        end
        
        
    end
end

