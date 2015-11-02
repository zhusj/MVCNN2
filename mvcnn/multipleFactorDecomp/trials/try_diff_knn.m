load('AVLetter_Learn_HOSVD_52-10-(1-2)');
load('results_NEW');
letter_bases = u{1};
num_train_reps = 2;
ltr_labels = (1:(length(letter_bases)/num_train_reps))';%[1:26,1:26]';
ltr_labels = repmat(ltr_labels,num_train_reps,1);

letter_act_labels = (1:26);%[1:26,1:26]';
letter_act_labels = repmat(letter_act_labels,10,1);
letter_act_labels  = reshape(letter_act_labels ,[],1);
for i = 1:5
    letter_est_results =  knnclassify(estimated_letter_vector,letter_bases,ltr_labels,i);
    letter_recog_rate = 1-(double(nnz(letter_est_results-letter_act_labels))/numel(letter_act_labels))
end
