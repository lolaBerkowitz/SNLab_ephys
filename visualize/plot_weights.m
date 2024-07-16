function plot_weights(ica,regionChan,basepath,nICs)

fig = figure;
fig.Color = [1 1 1];
plot(ica.topo,repmat([1:numel(ica.channels)]',1,nICs),'LineWidth',2);
set(gca,'YDir','reverse')
hold on
legend({num2str(ica.channels)})
ylim([1 size(ica.weights,2)])
if ~isempty(regionChan)
    for ii = 1:length(regionChan)
        idx = find(ica.channels == regionChan(ii));
        if ~isempty(idx)
            line([-300 300],[idx idx],'Color',[0.5 0.5 0.5],'LineStyle','--');
        end
    end
end
if ~isfolder('ICA')
    mkdir('ICA')
end
saveas(gcf,[basepath,filesep,'ICA\Weights.png']);

end