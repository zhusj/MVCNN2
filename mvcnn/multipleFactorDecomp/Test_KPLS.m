function [recog_rate,recog_rate_unit,unit_conflict] = Test_KPLS(pd,run_stamp,p_trained_data_file,p_split_data_file,dim,p_classifier,p_save_file)
% FILENAME: Test_KPLS.m
%train over 26x2x10, test for 26x1x10.
plot_results = 0;
if ~exist('p_save_file','var')
    p_save_file = 1;
end

classifier = 'RfC';%regression for classification
if ~isempty(strfind(p_classifier,'SVM'))  % Using SVM Classifier
    classifier = 'SVM';
end
%diary(sprintf('test_KPLS_log_%s.txt',run_stamp));
%echo on

p_use_wts=0;
p_use_sampling=0;

if ~exist('p_trained_data_file','var') || isempty(p_trained_data_file)
    return;
end
%MAP = KPLS;
load(p_trained_data_file);
%Data file has 'MAP', 'B_train','A_train','K_train','labels_train','new_labels'
MAP = MAP.setDimension(pd);
K_train = MAP.K;

%% learn SVM model for latent space
if strcmp(classifier,'SVM')
    svm_param_search = 0;
    %choose the best C_SVM parameters
    bestc=1; bestg=4;
    if svm_param_search
        fprintf('\n\nSVM grid parameter search:');
        pls_dim = size(MAP.T,2);
        est_g_value = mean(var(MAP.T))*pd; %mean variance of features*number of features;
        g_values = (0.5:0.05:1.5)*est_g_value;
        g_values = g_values(g_values>0);
        bestcv = 0;
        for log2c = -1:1:3,
          for log2g =-4:1:5
            cmd = ['-q -c ', num2str(2^log2c), ' -g ', num2str(2^log2g)];
            cv = get_cv_ac(train_labels, MAP.T, cmd, 3);
            if (cv >= bestcv),
              bestcv = cv; bestc = 2^log2c; bestg = 2^log2g;
            end
            sprintf('\n%g %g %g (best c=%g, g=%g, rate=%g)\n', log2c, log2g, cv, bestc, bestg, bestcv);
          end
        end
    end


    %SVM_model = ovrtrain(matrix_labels, MAP.T,sprintf('-c %g -g %g',bestc, bestg));
    SVM_model = ovrtrain(train_labels, MAP.T,sprintf('-c %g -g %g -q',bestc, bestg));
    %SVM_model2 = ovrtrain(labels_train, MAP.T,sprintf('-g 18'));
end

%%
if ~exist('p_split_data_file','var') || isempty(p_split_data_file)
    return;
end
load(p_split_data_file);

sz_test = size(test_GCF);
%prepare actual labels and data 
if numel(size(test_GCF)) > 2
    sz_test = sz_test(1:end-1);
    test_labels = get_unfolding_label(dim,sz_test);
    %prepare data tensor 
    A_test =  vectorize_CF(test_GCF);
else
    A_test = test_GCF;
end

if numel(size(train_GCF)) > 2
    A_train = vectorize_CF(train_GCF);
else
    A_train = train_GCF;
end

clear test_GCF
clear train_GCF




if strcmp(classifier,'SVM') % Using SVM Classifier
    effictive_formula = 'Kt*MAP.R2';
    effictive_formula2 = 'A_test*MAP.R';
else
    effictive_formula = 'Kt*MAP.M2';
    effictive_formula2 = 'A_test*MAP.M';
end
    


Kt = get_gram_kernel(A_test,A_train);
%Kt = At*A_train';
nt = size(A_test,1);
wt = ones(nt,1);
n = size(A_train,1);
w = ones(n,1);
wt= wt*w';
%w = eye(n)-w/n;
Kt = (Kt-wt/n*K_train)*(eye(n)-w*w'/n);%normalizing K using equation (13) in the paper
% 

test_results = eval(effictive_formula);%Kt*MAP.M2;
if strcmp(classifier,'SVM')
    [test_est_labels,acc]= ovrpredict(test_labels, test_results,SVM_model);
else
    [~,test_est_labels] = max(test_results,[],2);
end
diff = (test_est_labels-test_labels);
recog_rate = 1-double(nnz(diff))/numel(test_labels)

recog_rate2 = -1;
   
    

% %plotting(mapped_B_letter, letter_act_labels,1,ltr_fig);
% %plotting(mapped_B_style, style_act_labels,1,sty_fig);
% 
% style_est_results = knnclassify(mapped_B_style,style_bases,style_labels,5);


%% compute letter estimation error

%letter_est_labels = knnclassify(letter_est_results, letter_bases,labels_train,1,'cosine');
%letter_est_labels = mod(letter_est_labels-1,26)+1;
%letter_est_labels = classify_vector(letter_est_results, letter_bases,labels_train);



if numel(size(A_test))==2
    if size(A_test,1)>1
        A_test = bsxfun(@minus, A_test, mean(A_test));
    end
    test_results2 = eval(effictive_formula2);% At*MAP.M;
    if strcmp(classifier,'SVM')
        [test_est_labels2,acc]= ovrpredict(test_labels, test_results2,SVM_model);
    else
        [~,test_est_labels2] = max(test_results2,[],2);
    end

    diff = (test_est_labels2-test_labels);
    recog_rate2 = 1-double(nnz(diff))/numel(test_labels)

    if recog_rate2>recog_rate
        %temp1 = test_results;
        %temp2 = test_est_labels;
        temp3 = recog_rate;

        test_results = test_results2;
        test_est_labels = test_est_labels2;
        recog_rate = recog_rate2;
        effictive_formula = effictive_formula2;

        recog_rate2 = temp3;
    end
end

if length(test_labels)>10
    %figure,plot(test_labels);
    %hold on,plot(test_est_labels,'r.');


    recog_rate_unit = zeros(length(unique(test_labels)),1);
    for i = unique(test_labels)'
        indx = find(test_labels == i);
        d = diff(indx);
        recog_rate_unit(i) = 1-double(nnz(d))/numel(indx);
    end

    unit_conflict = zeros(length(unique(test_labels)));
    for i = 1:length(test_labels)
        unit_conflict(test_labels(i),test_est_labels(i)) = unit_conflict(test_labels(i),test_est_labels(i))+1;
    end
else
    recog_rate_unit = recog_rate;
    unit_conflict = [];
end
%plot_confusion_matrix(unit_conflict,


%% save data to file
if p_save_file
    save(sprintf('%s_%s',mfilename, run_stamp),'test_results','recog_rate2','recog_rate','recog_rate_unit','unit_conflict','effictive_formula','classifier');
end

% plotting(mapped_B_letter,letter_act_labels,1,ltr_fig,letter_diff==0);
% plotting(mapped_B_style,style_act_labels,1,sty_fig,style_diff==0);
%diary off
return 


    
