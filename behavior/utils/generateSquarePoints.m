function [x, y] = generateSquarePoints(corners, numPoints)
% GENERATESQUAREPOINTS Generate points along the sides of a square
%   Inputs:
%       corners - 4x2 matrix of [x,y] coordinates of square corners in order
%       numPoints - total number of points to generate (default 1000)
%   Outputs:
%       x, y - vectors of point coordinates

    % lazy variable input handling 
    if nargin < 2
        numPoints = 1000;
    end
    
    % Calculate the perimeter and side lengths
    sideLengths = zeros(4,1);
    for i = 1:3
        sideLengths(i) = norm(corners(i+1,:) - corners(i,:));
    end
    sideLengths(4) = norm(corners(1,:) - corners(4,:));
    perimeter = sum(sideLengths);
    
    % Calculate number of points per side (proportional to side length)
    pointsPerSide = round(numPoints * sideLengths / perimeter);
    
    % Adjust total points to match exactly 1000
    totalPoints = sum(pointsPerSide);
    while totalPoints ~= numPoints
        if totalPoints < numPoints
            [~, idx] = max(pointsPerSide);
            pointsPerSide(idx) = pointsPerSide(idx) + 1;
        else
            [~, idx] = min(pointsPerSide);
            pointsPerSide(idx) = pointsPerSide(idx) - 1;
        end
        totalPoints = sum(pointsPerSide);
    end
    
    % Generate points for each side
    x = [];
    y = [];
    
    for i = 1:4
        if i < 4
            startPt = corners(i,:);
            endPt = corners(i+1,:);
        else
            startPt = corners(4,:);
            endPt = corners(1,:);
        end
        
        % Parameter t goes from 0 to 1
        t = linspace(0, 1, pointsPerSide(i))';
        
        % Linear interpolation between start and end points
        sideX = startPt(1) + t * (endPt(1) - startPt(1));
        sideY = startPt(2) + t * (endPt(2) - startPt(2));
        
        % Avoid duplicate points at corners (except first/last point)
        if i > 1
            sideX = sideX(2:end);
            sideY = sideY(2:end);
        end
        
        x = [x; sideX];
        y = [y; sideY];
    end
end