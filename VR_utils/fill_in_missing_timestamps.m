
function interp_times = fill_in_missing_timestamps(timestamps,distances)
% given a set of incomplete timestamps (such as thosed used for syncing),
% and a set of distances between an independent timestamp source
% (distances), fill in the missing timestamps. 

% example 
% sync_ts = digitalIn.timestampsOn{1, 7}; % from intan
% distances = diff(vr_pos.experiment_ts); % from godot logs 
% interp_ts = fill_in_missing_timestamps(sync_ts,distances)




% Compute the number of intervals between timestamps
num_intervals = length(distances) + 1;

% Compute the time step for each interval
interval_times = [0; cumsum(distances)];
interval_steps = diff(interval_times);

% Compute the interpolated timestamps
interp_times = zeros(1, num_intervals);
interp_times(1) = timestamps(1);
for i = 2:num_intervals
    interp_times(i) = interp_times(i-1) + interval_steps(i-1);
end
end