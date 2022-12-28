function [outputArg1,outputArg2] = behaivor_summary(basepath)
%behaivor_summary provides a summary of general locomotion for behavior file
% found in basepath. 


basename = basenameFromBasepath(basepath);
session = loadSession(basepath,basename);

% for a given basepath, load behavior and session file. If no behavior file exists,
% return. 
try
    load(fullfile(basepath,[basename,'.animal.behavior.mat']));
catch
    error('Problem loading behavior file. Check basepath and try again') 
end

% identify epochs for behavior from session.behavioralTracking 
for i = 1:length(session.behavioralTracking)
    behav_ep(i,1) = session.epochs{1, session.behavioralTracking{1, i}.epoch}.startTime; 
    behav_ep(i,2) = session.epochs{1, session.behavioralTracking{1, i}.epoch}.stopTime; 
end

% Identify trials 
if ~isempty(behavior.trials)
    n_trials = size(behavior.trials,1);
else
    disp('No trials found. Computing for whole epoch')
    n_trials = NaN;
end

behavior_summary = table;

% For trials in epoch, get xy coordinates 

% path length while running 

% velocity while running 

% number and duration of stops 

% search area (overall occupancy) 

behavior


end

