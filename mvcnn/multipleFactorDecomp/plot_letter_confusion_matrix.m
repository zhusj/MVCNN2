function acc_grpd=plot_letter_confusion_matrix(p_letter_est_results,run_stamp,data_file)
    %LIST = [ 'A';'B';'C';'C';'C';'F';'G';'H';'I';'J';'A';'A';'M';'A';'O';'B';'O';'R';'S';'C';'O';'V';'W';'S';'Y';'C'];
    LIST = [ 'A';'B';'C';'C';'E';'F';'G';'H';'I';'G';'K';'L';'M';'N';'O';'B';'Q';'R';'S';'C';'Q';'V';'W';'S';'Y';'Z'];
    
    if ~exist('run_stamp','var')
        run_stamp = 'rep_3';
    end
    if ~isempty(p_letter_est_results)
        letter_est_results = p_letter_est_results;
    else
        if nargin ==2
            load(data_file,'letter_est_results');
        else
            acc_grpd = -1;
            return;
        end
    end
    numLetter = length(LIST);
    letter_confusion = zeros(numLetter,numLetter);
    for i = 1:numLetter
        for j = 1:numLetter
            letter_confusion(i,j) = sum(letter_est_results(i,:)==j);
        end
    end
    acc = sum(diag(letter_confusion))/sum(sum(letter_confusion))
    %set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
    figure, imagesc(letter_confusion);
    %colormap(gray);
    %axis(['A' 'Z' 'A' 'Z']);
    axis_labels = {['A':'Z']'};
    title('Letter Not Grouped');
    set(gca,'XTick',1:length(axis_labels));
    set(gca,'XTickLabel',axis_labels);
    set(gca,'YTick',1:length(axis_labels));
    set(gca,'YTickLabel',axis_labels);
    grid on
    
%     numGroups = length(unique(LIST));
%     letter_confusion_grpd = zeros(numGroups,numGroups);
    letter_confusion_grpd = zeros(numLetter,numLetter);
    for i = 1:numLetter
        for j = 1:numLetter
            mapped_i = LIST(i) - 'A' +1;
            mapped_j = LIST(j) - 'A' +1;
            letter_confusion_grpd(mapped_i,mapped_j) = letter_confusion_grpd(mapped_i,mapped_j) + sum(letter_est_results(i,:)==j);
            %letter_confusion(i,j) = sum(letter_est_results(i,:)==j);
        end
    end
    acc_grpd = sum(diag(letter_confusion_grpd))/sum(sum(letter_confusion_grpd))
    %set(gca,'XTickLabel',{'-pi','-pi/2','0','pi/2','pi'})
    figure, imagesc(letter_confusion_grpd);
    %colormap(gray);
    %axis(['A' 'Z' 'A' 'Z']);
    title('Letter Grouped');
    set(gca,'XTick',1:length(axis_labels));
    set(gca,'XTickLabel',axis_labels);
    set(gca,'YTick',1:length(axis_labels));
    set(gca,'YTickLabel',axis_labels);
    grid on
    save(sprintf('letter_confusion_%s',run_stamp),'letter_confusion','letter_confusion_grpd');
end
