function fig = plot_path(behavior,session)
% plots path in behavior.position by trial.
% LB 12/2022

if ~isfield(behavior,'trials')
    error('No trials field found in behavior file')
end

% behavior epochs 
behave_ep = [];
for ep = 1:length(session.behavioralTracking)
    behave_ep(ep,:) = [session.epochs{1,session.behavioralTracking{1,ep}.epoch}.startTime, ...
        session.epochs{1,session.behavioralTracking{1,ep}.epoch}.stopTime];
end
behave_ep = IntervalArray(behave_ep);
trial_ep = IntervalArray(behavior.trials);

fig = figure;
for i = 1:behave_ep.n_intervals
    
    cur_ep = behave_ep(i) & trial_ep;
    idx = InIntervals(behavior.timestamps, cur_ep.intervals);
 
    x = behavior.position.x(idx);
    y = behavior.position.y(idx);
    [occ,~] = behavior_funcs.occ_map(x,y,2,behavior.sr);
    axs(1) = subplot(2,behave_ep.n_intervals,i,'align');
    b = imagesc(flipud(occ));
    axis image
    axis equal
    axis off
    set(b,'AlphaData',~isnan(flipud(occ)))
        title(behavior.trialsID{i})

    axs(2) = subplot(2,behave_ep.n_intervals,i+behave_ep.n_intervals) ;
    plot(x,y,'k')
    axis image
    axis equal
    axis off
end

end