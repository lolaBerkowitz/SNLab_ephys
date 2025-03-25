function [xOuter, yOuter, xInner, yInner] = generateNestedSquares(basepath, ep,varargin)
% GENERATENESTEDSQUARES Generate points on outer and inner square sides
%   Inputs:
%       basepath: location of behavioral data kept in CellExplorer format. 
%       ep: integer indicated behavior epoch that contains maze_coords in session.behavioral_tracking. 
%       numPoints - total number of points per square (default 1000)
%       scaleFactor - scaling factor for inner square (default 0.8)
%   Outputs:
%       xOuter, yOuter - vectors of outer square point coordinates
%       xInner, yInner - vectors of inner square point coordinates

    p = inputParser;
    addParameter(p,'numPoints',1000);
    addParameter(p,'scaleFactor',.8);

    % type=parse(p,Parameters,'');
    p.parse(varargin{:});

    numPoints=p.Results.numPoints;
    scaleFactor=p.Results.scaleFactor;

    % load session file 
    basename = basenameFromBasepath(basepath); 
    session = loadSession(basepath,basename); 
    
    % load maze coordinates for the input epoch
    maze_coords = session.behavioralTracking{1,ep}.maze_coords; 
    corners = [maze_coords.x_scaled, maze_coords.y_scaled];
    
    % Generate points for outer square
    [xOuter, yOuter] = generateSquarePoints(corners, numPoints);
    
    % Calculate center of outer square
    center = mean(corners, 1);
    
    % Create inner square (scaled down and centered)
    innerCorners = zeros(size(corners));
    for i = 1:4
        % Vector from center to corner
        vec = corners(i,:) - center;
        % Scale the vector
        scaledVec = vec * scaleFactor;
        % New corner position
        innerCorners(i,:) = center + scaledVec;
    end
    
    % Generate points for inner square
    [xInner, yInner] = generateSquarePoints(innerCorners, numPoints);
end