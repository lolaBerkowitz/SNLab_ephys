function update_metadata_pixel_distance(metadata_path)

df = readtable(metadata_path,'Delimiter','comma'); 
for i = 1:length(df.basepath)
    
    if ~isnan(df.pixel_distance(i))
        continue
    end
    
    if ismember(df.vidname{i},{'na','MISSING'})
        df.pixel_distance(i) = nan;
    else
        try
            video_path = fullfile(df.basepath{i},[df.vidname{i},'.avi']);
            pixel_reference = df.pixel_dist_reference(i);
            
            pixel_distance = tracking.maze_distance_gui(video_path,pixel_reference);
    
            df.pixel_distance(i) = pixel_distance; 
            
            
            writetable(df,metadata_path)
        catch 
            video_path = fullfile(df.basepath{i},[df.vidname{i},'.mp4']);
            pixel_reference = df.pixel_dist_reference(i);
            
            pixel_distance = tracking.maze_distance_gui(video_path,pixel_reference);
    
            df.pixel_distance(i) = pixel_distance; 
            
            
            writetable(df,metadata_path)
        end
    end
      
end

end