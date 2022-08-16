function [ angle ] = xy_angle(x,y)
%XYangle Obtains angle of heading from XY path*
%  Input
%       x: matrix containing x(1), x(2)
%       y: matrix containing y(1), y(2)
%  Output
%       angle: single column vector of angle of heading
% 
% *Main calculation: angle= atan2d(y2-y1,x2-x1) + 360*((y2-y1)<0);
% 
% Ryan Harvey 6/2/17

    angle = deg2rad(atan2d(y(:,2)-y(:,1),x(:,2)-x(:,1)) + 360*((y(:,2)-y(:,1))<0));
end


