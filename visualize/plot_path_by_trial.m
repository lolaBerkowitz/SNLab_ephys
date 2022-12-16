function fig = plot_path_by_trial(behavior)
% plots path in behavior.position by trial.
% LB 12/2022

if ~isfield(behavior,'trials')
    error('No trials field found in behavior file')
end

fig = figure;
for i = 1:length(behavior.trialsID)
    idx = behavior.timestamps > behavior.trials(i,1)...
        & behavior.timestamps < behavior.trials(i,2);
    subplot(1,length(behavior.trialsID),i)
    plot(behavior.position.x(idx),behavior.position.y(idx),'k')
    title(behavior.trialsID{i})
    axis image
    axis equal
    axis off
end

end