function plot_weights_individual(ica,basepath,nICs)
fig=figure('DefaultAxesFontSize',8,'defaultAxesFontName','Serif','defaultTextFontName','Serif');
fig.Color = [1,1,1];
[fig_width_in, fig_height_in] = set_size('paper', .75, [1,1]);
set(fig,'Position',[1 1 fig_width_in fig_height_in])

for i = 1:nICs
    subplot(1,nICs,i)
    plot(ica.topo(:,i),repmat([1:numel(ica.channels)]'',1,1),'LineWidth',2,'Color','k');
    hold on;
    plot(ones(numel(ica.channels)),repmat([1:numel(ica.channels)]',1,1),'Color','r')
    if i == 1
        set(gca,'YDir','reverse','ytick',repmat([1:numel(ica.channels)]',1,1),'yticklabel',ica.region)
    else
        set(gca,'YDir','reverse')
    end
    hold on
    title({'variance: ', num2str(round(ica.meanvar(i)))})
    ylim([1 size(ica.weights,2)])
end

if ~isfolder('ICA')
    mkdir('ICA')
end
% saveas(gcf,[basepath,filesep,'ICA\individual_weights.png']);
end