df = readtable('Y:\laura_berkowitz\behavior_metadata.csv','Delimiter','comma'); 

for i = 1:length(df.basepath)
    
    if ~isnan(df.pixel_distance(i))
        continue
    end
    
    if ismember(df.vidname{i},'na')
        df.pixel_distance(i) = nan;
    else
        video_path = fullfile(df.basepath{i},[df.vidname{i},'.avi']);
        maze_size = df.maze_width_cm(i);
        
        pixel_distance = tracking.maze_distance_gui(video_path,maze_size,300);

        df.pixel_distance(i) = pixel_distance; 
        
        writetable(df,'Y:\laura_berkowitz\behavior_metadata.csv')
    end
    
    
end