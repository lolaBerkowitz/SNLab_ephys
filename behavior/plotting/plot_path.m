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
fig.Color = [1 1 1];
for i = 1:behave_ep.n_intervals
    maze_coords = session.behavioralTracking{1, i}.maze_coords;
    obj_idx = contains(maze_coords.object,{'A','B'}) & contains(maze_coords.position,{'center'});
    boundary_idx = contains(maze_coords.object,{'corner'});

    cur_ep = behave_ep(i) & trial_ep;
    idx = InIntervals(behavior.timestamps, cur_ep.intervals);
 
    x = behavior.position.x(idx);
    y = behavior.position.y(idx);
    [occ,~] = behavior_funcs.occ_map(x,y,2,60);
    axs(1) = subplot(2,behave_ep.n_intervals,i,'align');
    b = imagesc(flipud(occ));
    axis image
    axis equal
    axis off
    set(b,'AlphaData',~isnan(flipud(occ)))
        title(behavior.trialsID{i})
    hcb = colorbar;
    set(get(hcb,'Title'),{'String','Rotation','Position'},{'Seconds',90,[35 50]})
    
    
    axs(2) = subplot(2,behave_ep.n_intervals,i+behave_ep.n_intervals) ;
    plot(x,y,'k'); hold on;
    scatter(maze_coords.x_scaled(obj_idx),maze_coords.y_scaled(obj_idx),'*r');
    scatter(maze_coords.x_scaled(boundary_idx),maze_coords.y_scaled(boundary_idx),'*b');
    hold off
    legend({'path','objects','maze boundaries'},'Position',[0.5,0.01,0.15,0.075])
    legend('boxoff') 
    axis image
    axis equal
    axis off
    
end

end