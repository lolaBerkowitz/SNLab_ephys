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
    [occ,~] = behavior_funcs.occ_map(behavior.position.x(idx),behavior.position.y(idx),3,behavior.sr);
    axs(1) = subplot(2,length(behavior.trialsID),i,'align');
    b = imagesc(flipud(occ));
    axis image
    axis off
    set(b,'AlphaData',~isnan(flipud(occ)))
        title(behavior.trialsID{i})

    axs(2) = subplot(2,length(behavior.trialsID),i+length(behavior.trialsID)) ;
    plot(behavior.position.x(idx),behavior.position.y(idx),'k')
    axis image
    axis off
end

end