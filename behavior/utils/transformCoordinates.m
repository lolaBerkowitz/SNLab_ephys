function [ x_cm,y_cm] = transformCoordinates(behavior) 
%transformCoordinates converts xy coordinates using max and min input
%values by translating and then scaling (by maze diameter). 
%   input: 
%       - dia: maze diameter (assumes symetric environment) in cm
%       - x_max: maze x_max in pixels
%       - x_min: maze x_min in pixels
%       - y_max: maze y_max in pixels
%       - y_min: maze y_min in pixels
%       - x: x coordinates in pixels
%       - y: y coordinates in pixels
%   output: 
%       -x and y: coordinates centered at 0,0 and scaled to centimeters. 
    
    % load session
    x = behavior.position.x;
    y = behavior.position.y;
    
    % for a given epoch 
    maze_origin_x = median(min(maze_coords.x_cm(ismember(maze_coords.object,'corner'))):...
        max(maze_coords.x_cm(ismember(maze_coords.object,'corner'))));

    maze_origin_y = median(min(maze_coords.y_cm(ismember(maze_coords.object,'corner'))):...
        max(maze_coords.y_cm(ismember(maze_coords.object,'corner'))));
    %Find Origin
    x_Origin = median(min(x):max(x)); 
    y_Origin = median(min(y):max(y));
      
    %% Convert tracking data to cm
    x_cm =((x-(x_Origin))*(dia /(x_max - x_min)));
    y_cm =((y-(y_Origin))*(dia /(y_max - y_min)));
end

