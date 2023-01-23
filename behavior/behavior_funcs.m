classdef behavior_funcs
    %behavior consists of functions used to assess locomotion of rodent
    % behavior
    %
    % Adapted from OF from ephys_tools by LB 12/2022
    
    
    
    properties
        Property1
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
            length = sum(distance_vector(velocity >= run_threshold,1));%Total Path length for points greater than 3cm/s
            
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
        
        function stop_measures = stops(x,y,ts,velocity,fr,epoch,varargin)
            % Finds when animal stop ( < 3cm/s velocity for at least 1
            % second) and computes stop features.
            % inputs:
            %       x: position vector for x of length n
            %       y: position vector for y of length n
            %       ts: timestamps of length n
            %       fr: frame rate
            %       velocity: instantaneous velocity of length n-1
            %       epoch: length of minimum stop epoch in frames
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
            
            parse(p,varargin{:});
            run_threshold = p.Results.run_threshold;
            
            [stop_measures.stopIdx,startStop,endStop] = behavior_funcs.stop_epochs(velocity,epoch,'run_threshold',run_threshold);
            
            %     This finds coordinates associated with stop
            for ii=1:length(startStop)
                motionless{ii} = [x(startStop(ii):endStop(ii)),y(startStop(ii):endStop(ii))];
                timeMotionless{ii} = size(motionless{ii},1)/fr;
                tsStop{ii} = ts(startStop(ii):endStop(ii));
            end
            
            stop_measures.stops = motionless;
            stop_measures.timeStopped = timeMotionless;
            stop_measures.tsStop = tsStop;
            
            %Create Number of Stops
            stop_measures.NumStops = size(stop_measures.stops,2);
            
            %Find center of mass for stops.
            for ii=1:size(stop_measures.stops,2)
                center(ii) = center_of_mass(stop_measures.stops{1,ii}(:,1),stop_measures.stops{1,ii}(:,2));
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
            
            HD = wrapTo360(fixNLXangle(rad2deg(atan2(head(:,2)-nose(:,2),...
                head(:,1)-nose(:,1))),round(0.1667*30)));
            angVel = insta_angvel(HD',fr);
            
            
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
            
            xedge=linspace(min(x),max(x),round(range(x)/binsize));
            yedge=linspace(min(y),max(y),round(range(y)/binsize));
            
            %Bin coordinates for heat map and apply guassian filter to smooth
            [map] = histcounts2(x,y,xedge,yedge);
            map=flipud(imrotate(map,90));
            map=map/fr;
            
            % smooths over 1.5 cm
            occ = imgaussfilt(map, 1);
            occ(map==0) = nan;
            map(map==0) = nan;
        end
        
        function [out] = thigmotaxis(x,y,fr,diameter,varargin)
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
            
            
            sa = sum(sum(map>0))/sum(sum(~isnan(map))); %Calculate the proportion of bins occupied by animal.
            
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
        
        %% Object exploration functions
        
        function [bound_x, bound_y] = object_boudary(object_center,object_edge)
            % Given the center and edge of an object, compute and return xy
            % coordinates representing the radius around an object
            % input:
            % object_center: [x,y] from maze_coords.csv
            % object_edge: [x, y]
            % center of object
            x0 = object_center(:,1);%%origin
            y0 = object_center(:,2);%%origin
            
            % radius (difference from center to edge
            r = abs(sqrt((object_center(1,1)-object_edge(1,1))^2+(object_center(1,2)-object_edge(1,2)).^2));
            
            % create circle of 25 points around object
            bound_x = (cos(linspace(-pi,pi,25)) * r) + x0;
            bound_y = (sin(linspace(-pi,pi,25)) * r) + y0;
            
        end
        
        function metrics = object_explore_vectors(object_center,object_edge,behavior,varargin)
            % calculate object interaction for all xy coordinates and
            % return as structure.  LB 2022
            
            p = inputParser;
            p.addParameter('perimeter',3,@isnumeric) % perimeter around object in units as coordinates (default cm)
            p.addParameter('run_thres',5,@isnumeric)
            
            p.parse(varargin{:})
            perimeter = p.Results.perimeter;
            run_thres = p.Results.run_thres;
            
            % grab xy coordinates
            x = behavior.position.x;
            y = behavior.position.y;
            
            % compute velocity
            [velocity, ~,~] = linear_motion(x,y,behavior.sr,1);
            
            % get object boundary
            [bound_x, bound_y] = object_boudary(object_center,object_edge,perimeter);
            
            % logical index for all coordinates inside object boundary
            in_home = inpolygon(x,y,bound_x',bound_y');
            
            % time moving slow for at least 1 second
            stopIdx = contiguousframes(velocity < run_thres,behavior.sr);
            
            % coordinates for frames inside home base for at least 0.5 seconds
            out_home = contiguousframes(in_home,behavior.sr/2); %has to be inside of hb for at least 1 sec to count as entryd
            
            % This finds stops that occur in the home base boundary
            [startStop,~,~] = find_groups(stopIdx);
            
            % get n entries into object area
            [~,~,metrics.entries] = find_groups(out_home);
            
            % total time near object
            metrics.obj_Occ = nansum(out_home)/fr;
            
            % average velocity near object
            metrics.obj_Vel = nanmean(velocity(out_home(1:end-1,1),1)); %remove last tempIn idx to accomodate velocity length
            
            % time moving slow near object in seconds
            metrics.obj_slow = nansum(out_home(1:end-1) & stopIdx)/fr; % time being slow in homebase
            
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
        
%         %% locomotion over time 
%         function locomotion_over_epochs(x,y,[startIdx,endIdx])
%             % examine path length, search area, velocity, and stops as a
%             % function of epoch (could be trials, behavior epochs, or time
%             % bins
%             disp('NOT FUNCTIONTIONAL LB 12/2022')
%         end
%         
        %%  behavior utils
        
        
    end
end

