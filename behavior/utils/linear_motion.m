function [velocity, acceleration, distance_vector] = linear_motion(x,y,t,fr)
% computes linear velocity and acceleration from position coordinates
% inputs: 
%   x: position vector of x-coordinates of length n
%   y: position vector of y-coordinates of length n
%   t: position vectpr of timestamps of length n
%   fr : video frame rate in Hz 
% output:
%   velocity: n-1 vector of velocity
%   acceleration: second derivative of position data using matlab gradient
%   function. 
%  LB 2020
%squared distance between consecutive points
sqrXDiff = (diff(x)).^2;
sqrYDiff = (diff(y)).^2;

%distance formula
distance_vector = sqrt(sqrXDiff + sqrYDiff);
velocity = LinearVelocity([t, x, y]);
velocity = velocity(:,2);
acceleration = gradient(velocity,1/fr); %instanteous acceleration

end
