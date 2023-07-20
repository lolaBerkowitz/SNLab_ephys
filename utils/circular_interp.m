function time_phase =  circular_interp(ts,theta,times)
% interpolates spike times to angle data using sin/cosine 
%input:
%  ts: timestamp vector (length n)
%  theta: angle in degress (length n)
%  spks:  times 
%output: 
%  spk_phase: angles in degrees when spike occured
% get sin and cos then interp spike times
time_phase = atan2(interp1(ts,sin(deg2rad(theta)),times,'linear'),interp1(ts,cos(deg2rad(theta)),times,'linear'));
% wrap to 360
time_phase = wrapTo360(rad2deg(time_phase));

end