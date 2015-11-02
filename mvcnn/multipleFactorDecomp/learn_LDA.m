function storeFile = learn_LDA(filename,trainRepeats ,plot_flg)
% FILE NAME: HOFDA.m (Higher Order Fisher Diecreminating Analysis)
% tensor analysis for gcf not gscf

    diary 'run_log_learn_HOFDA.txt'
    echo on

    if nargin<2
        plot_flg = 1;
    end

    if nargin<1
        filename = 'AVLetter_D1-Letters_D2-reptns(1-2)_D3-Subjects_D4-5-SequenceFeatures';
    end
    load(filename);
    pack
    
    %prepare labels
    letter_labels = (1:size(GCF,1))';%zeros(1,nSeqs);
    %letter_labels = (1:(size(GCF,1)*size(GCF,2)))';%zeros(1,nSeqs);
    letter_labels = repmat(letter_labels,size(GCF,2)*size(GCF,3),1);
    %letter_labels = repmat(letter_labels,size(GCF,3),1);
    style_labels = 1:size(GCF,3);%zeros(1,nSeqs);
    style_labels = repmat(style_labels,size(GCF,1)*size(GCF,2),1);
    style_labels = reshape(style_labels,[],1);
    
    %prepare the data matrix
    s = size(GCF);
    B = reshape(GCF,prod(s(1:3)),[]);
    
    clear GCF;
    %% start LDA processing
    [letter_bases,MAP_L] = lda(B,letter_labels,length(unique(letter_labels))-1);
    [style_bases,MAP_P] = lda(B,style_labels,length(unique(style_labels))-1);
    letterMapping = MAP_L.M;
    styleMapping = MAP_P.M;
    letter_mean = MAP_L.mean;
    style_mean = MAP_P.mean;
    
    letter_labels = mod(letter_labels-1,26)+1;
    
    storeFile = sprintf('AVLetter_Learn_LDA_52-10-(%s)',cell2string(trainRepeats));
    save(storeFile,'letter_bases','letterMapping','letter_labels','style_bases','styleMapping','style_labels','letter_mean','style_mean');
    %% Plotting
    if plot_flg
%         plotting(B,letter_labels);
%         plotting(B,style_labels);
        plotting(letter_bases,letter_labels);
        plotting(style_bases,style_labels);
    end
    echo off
    diary off
end
        
