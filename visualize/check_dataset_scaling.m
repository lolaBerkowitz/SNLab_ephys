
function check_dataset_scaling(dataset_path)

% Compile sessions
df = compile_sessions(dataset_path);

K = length(df.basepath);
ncols = 10;
nrows = ceil(K / ncols);

% --- First pass: get global axis limits ---
all_x = [];
all_y = [];

for k = 1:K
    basepath = df.basepath{k};
    basename = basenameFromBasepath(basepath);

    load(fullfile(basepath, [basename, '.animal.behavior.mat']), 'behavior');

    all_x = [all_x; behavior.position.x(:)];
    all_y = [all_y; behavior.position.y(:)];
end

% Remove NaNs just in case
all_x = all_x(~isnan(all_x));
all_y = all_y(~isnan(all_y));

xlims = [min(all_x), max(all_x)];
ylims = [min(all_y), max(all_y)];

% --- Plot ---
figure;

for k = 1:K
    subplot(nrows, ncols, k)

    basepath = df.basepath{k};
    basename = basenameFromBasepath(basepath);

    load(fullfile(basepath, [basename, '.animal.behavior.mat']), 'behavior');

    plot(behavior.position.x, behavior.position.y, 'k-.','LineWidth',.5)
    rectangle('Position',[xlims(1), ylims(1), diff(xlims), diff(ylims)], ...
        'EdgeColor','k','LineWidth',0.5)
    % enforce identical axes
    xlim(xlims)
    ylim(ylims)

    axis equal
    axis tight

    title(strrep(basename, '_', '\_'), 'fontsize', 6)
end


end
