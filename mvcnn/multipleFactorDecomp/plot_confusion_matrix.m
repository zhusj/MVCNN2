function plot_confusion_matrix(confusion_matrix,labels,run_stamp,save_plot)
    % load('configuration.mat')
    % plot_confusion_matrix(unit_conflict,database.dimensions(1).values)
    if ~exist('run_stamp','var')
        c = clock;
        run_stamp = mat2str(c(end-3:end));
    end
    if ~exist('save_plot','var')
        save_plot = 0;
    end
    
    %set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
    h = figure; imagesc(confusion_matrix);
    %colormap(gray);
    
    axis_labels = labels;
    title(sprintf('Confusion Matrix (%s)',run_stamp));
    set(gca,'XTick',1:length(axis_labels));
    set(gca,'XTickLabel',axis_labels);
    set(gca,'YTick',1:length(axis_labels));
    set(gca,'YTickLabel',axis_labels);
    grid on
    
    if save_plot
        print(h, '-dpng', sprintf('confusion_plot_%s',run_stamp));
    end
end
