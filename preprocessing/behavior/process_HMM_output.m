% df = readtable("U:\users\ryanh\projects\hpc_ctx\sessions_v2.csv",'Delimiter',',');
% df = df(contains(df.use,'TRUE'),:);
df = table()
df.basepath = {pwd};

check_original_states = true;
% visual of current edges (states) for linearized pos
% gui labeled as this:
%    ___ ___
%   | 1 | 4 |
%   |2  |0  |5
%   |__ | __|
%     3   6
% becomes this when linearized:
%    ___ ___
%   | 3 | 4 |
%   |5  |0  |6
%   |__ | __|
%     1   2
% This script converts to this:
%    ___ ___
%   |   |   |
%   |   |0  |
%   |__ | __|
%   2       1
%
% This script will also bring down new state 2 to same level as state 1

for i = 1:length(df.basepath)
    disp(df.basepath{i})
    run(df.basepath{i}, check_original_states)
end

function run(basepath, check_original_states)

basename = basenameFromBasepath(basepath);
if ~exist(fullfile(basepath, [basename, '.animal.behavior.mat']), 'file')
    return
end
load(fullfile(basepath, [basename, '.animal.behavior.mat']), 'behavior')

if ~isfield(behavior, 'states')
    return
end

epochs = load_epoch('basepath', basepath);
if all(cellfun(@isempty, epochs.environment))
    disp('add environments')
    gui_session(basepath)
    epochs = load_epoch('basepath', basepath);
end
% select just wmaze
epochs = epochs(contains(epochs.environment, {'tmaze', 'Mwheel', 'Tmaze'}), :);

% locate tmaze timestamps
epoch_ts_idx = InIntervals(behavior.timestamps, [epochs.startTime, epochs.stopTime])';


states = unique(behavior.states(epoch_ts_idx));
states = states(~isnan(states));
if length(states) == 3 || isempty(states)
    return
end

% figure to verify states
if check_original_states
    idx = contains(epochs.environment, {'tmaze', 'Mwheel', 'Tmaze'});
    for ind = find(idx)'
        [idx, ~, ~] = InIntervals(behavior.timestamps, [epochs.startTime(ind), epochs.stopTime(ind)]);
        figure;
        states = unique(behavior.states(idx));
        states = states(~isnan(states));
        for s = states
            scatter(behavior.position.x(idx' & behavior.states == s), ...
                behavior.position.y(idx' & behavior.states == s))
            hold on
        end
    end
    prompt = "Do the initial states look good? y/n: ";
    txt = input(prompt, "s");
    if contains(txt, 'n')
        disp('try turning HMM off and re-run')
        return
    end
end

% drop states 2 and 4 down to the level of the other two arms (1&3) so
% right and left will be contingous
idx = (behavior.states == 4 | behavior.states == 6 | behavior.states == 2) & epoch_ts_idx;

behavior.position.linearized(idx) = behavior.position.linearized(idx) - ...
    min(behavior.position.linearized(idx)) + ...
    max(behavior.position.linearized(behavior.states == 0));

% relabel states
states = behavior.states;
idx = (behavior.states == 4 | behavior.states == 6 | behavior.states == 2) & epoch_ts_idx;
states(idx) = 2;

idx = (behavior.states == 3 | behavior.states == 5 | behavior.states == 1) & epoch_ts_idx;
states(idx) = 1;
behavior.states = states;

% update statenames
behavior.stateNames = {'middle', 'right', 'left'};

idx = contains(epochs.environment, {'tmaze', 'Mwheel', 'Tmaze'});
tmaze_epochs = [epochs.startTime(idx), epochs.stopTime(idx)];
[idx, ~, ~] = InIntervals(behavior.timestamps, tmaze_epochs);
if sum(idx) == 0
    return
end

% figure to verify states
colors = {[0, 0.4470, 0.7410], [0.8500, 0.3250, 0.0980], [0.9290, 0.6940, 0.1250]};
figure;
subplot(1, 2, 1)
states = unique(behavior.states(epoch_ts_idx));
states = states(~isnan(states));
for s = states
    scatter(behavior.position.x(epoch_ts_idx & behavior.states == s), ...
        behavior.position.y(epoch_ts_idx & behavior.states == s), [], colors{s+1})
    hold on
end
title('whole first epoch')

subplot(1, 2, 2)
states = unique(behavior.states(epoch_ts_idx));
states = states(~isnan(states));
for s = states
    scatter(behavior.timestamps(epoch_ts_idx & behavior.states == s), ...
        behavior.position.linearized(epoch_ts_idx & behavior.states == s), [], colors{s+1})
    hold on
end
xlim([tmaze_epochs(1, 1), tmaze_epochs(end, end)])
% title('first 5 minutes linear')


prompt = "Does this look good? y/n: ";
txt = input(prompt, "s");
if contains(txt, 'y')
    save(fullfile(basepath, [basename, '.animal.behavior.mat']), 'behavior')
end

end