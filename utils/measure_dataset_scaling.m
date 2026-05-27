
function measure_dataset_scaling(dataset_path)

df = compile_sessions(dataset_path);

K = length(df.basepath);

aspect_ratios = nan(K,1);
x_ranges = nan(K,1);
y_ranges = nan(K,1);

for k = 1:K
    basepath = df.basepath{k};
    basename = basenameFromBasepath(basepath);

    load(fullfile(basepath, [basename, '.animal.behavior.mat']), 'behavior');

    x = behavior.position.x;
    y = behavior.position.y;

    % remove NaNs
    valid = ~isnan(x) & ~isnan(y);
    x = x(valid);
    y = y(valid);

    x_range = max(x) - min(x);
    y_range = max(y) - min(y);

    x_ranges(k) = x_range;
    y_ranges(k) = y_range;

    aspect_ratios(k) = x_range / y_range;
end

figure;
subplot(3,1,1)
plot(x_ranges,'o-'); ylabel('x range (cm)')

subplot(3,1,2)
plot(y_ranges,'o-'); ylabel('y range (cm)')

subplot(3,1,3)
plot(aspect_ratios,'o-'); ylabel('x/y aspect ratio'); xlabel('session')
yline(median(aspect_ratios),'r--')

end