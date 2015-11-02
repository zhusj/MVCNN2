function [recog_rate,recog_rate_unit] = Test_KPLS_SVM(run_stamp,p_trained_data_file,p_split_data_file,dim)
% FILENAME: Test_KPLS.m
%train over 26x2x10, test for 26x1x10.


%diary(sprintf('test_KPLS_log_%s.txt',run_stamp));
%echo on

p_use_wts=0;
p_use_sampling=0;

load(p_split_data_file);

sz_test = size(test_GCF);
sz_test = sz_test(1:end-1);

%prepare actual labels
test_act_labels = get_unfolding_label(dim,sz_test);

%prepare data tensor
B = reshape(test_GCF,prod(sz_test),[]);
%X = bsxfun(@minus, X, mapping.mean);
clear test_GCF;


if ~exist('p_trained_data_file','var') || isempty(p_trained_data_file)
    return;
end
load(p_trained_data_file);
%Data file has 'MAP', 'B_train','A_train','K_train','labels_train','new_labels'

%embed test data on new reduced space
% mapped_B_letter = bsxfun(@minus, B, letter_mean)*letterMapping;
% mapped_B_style = bsxfun(@minus, B, style_mean)*styleMapping;

At =  B;%get_gram_kernel(B,B_train);
% %At = bsxfun(@minus, At, mean(At));
%K_train = exp(-0.5*pdist2(B_train,B_train));
%
Kt = get_gram_kernel(At,A_train);
%Kt = At*A_train';
nt = size(B,1);
wt = ones(nt,1);
n = size(A_train,1);
w = ones(n,1);
wt= wt*w';
%w = eye(n)-w/n;
Kt = (Kt-wt/n*K_train)*(eye(n)-w*w'/n);%normalizing K using equation (13) in the paper
% 

test_results = Kt*MAP.R2;



% %plotting(mapped_B_letter, letter_act_labels,1,ltr_fig);
% %plotting(mapped_B_style, style_act_labels,1,sty_fig);
% 
% style_est_results = knnclassify(mapped_B_style,style_bases,style_labels,5);


%% compute letter estimation error

%letter_est_labels = knnclassify(letter_est_results, letter_bases,labels_train,1,'cosine');
%letter_est_labels = mod(letter_est_labels-1,26)+1;
%letter_est_labels = classify_vector(letter_est_results, letter_bases,labels_train);

%[test_est_labels,acc]= ovrpredict(get_matrix_label(test_act_labels,max(test_act_labels)), test_results,SVM_model);
%[~,test_est_labels] = max(test_results,[],2);
[test_est_labels,acc]= ovrpredict(test_act_labels, test_results,SVM_model1);
acc
diff = (test_est_labels-test_act_labels);
recog_rate = 1-(double(nnz(diff))/numel(test_act_labels))

recog_rate2 = -1;
effictive_formula = 'Kt*MAP.R2';
% if exist('SVM_model2','var')
%     [est_labels2,acc2] = ovrpredict(test_act_labels, test_results,SVM_model2);
% else
%     [est_labels2,acc2] = ovrpredict(test_act_labels, test_results2,SVM_model1);
% end
if numel(size(At))==2
    At = bsxfun(@minus, At, mean(At));
    test_results2 = At*MAP.R;
    [est_labels2,acc2] = ovrpredict(test_act_labels, test_results2,SVM_model1);
    acc2
    %[~,est_labels2] = max(test_results2,[],2);
    diff = (est_labels2-test_act_labels);
    recog_rate2 = 1-double(nnz(diff))/numel(test_act_labels)
end

if recog_rate2>recog_rate
    effictive_formula = 'At*MAP.R';
    %temp1 = test_results;
    %temp2 = test_est_labels;
    temp3 = recog_rate;
    
    test_results = test_results2;
    test_est_labels = est_labels2;
    recog_rate = recog_rate2;
    
    recog_rate2 = temp3;
end


figure,plot(test_act_labels);
hold on,plot(test_est_labels,'r.');

recog_rate_unit = zeros(length(unique(test_act_labels)),1);
for i = unique(test_act_labels)'
    indx = find(test_act_labels == i);
    d = diff(indx);
    recog_rate_unit(i) = 1-(double(nnz(d))/numel(indx));
end

unit_conflict = zeros(length(unique(test_act_labels)));
for i = 1:length(test_act_labels)
    unit_conflict(test_act_labels(i),test_est_labels(i)) = unit_conflict(test_act_labels(i),test_est_labels(i))+1;
end



%% save data to file
save(sprintf('TEST_KPLS_SVM_%s',run_stamp),'test_results','recog_rate2','recog_rate','recog_rate_unit','unit_conflict','effictive_formula');

% plotting(mapped_B_letter,letter_act_labels,1,ltr_fig,letter_diff==0);
% plotting(mapped_B_style,style_act_labels,1,sty_fig,style_diff==0);
%diary off
return 


    
