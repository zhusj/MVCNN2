function pls_dim = get_PLS_dim_CV(A_train,labels_train,split_labels)
%function for  computing best PLS_dim using cross-validation
%global database

matrix_labels = get_matrix_label(labels_train);
n = size(A_train,1);
v_thresh = 1e-5;
m_thresh = 1e-3;
nfolds = 5;
temp = randperm(n);
rand_ind = temp;%(randperm(n));%double step random generation for ensuring complete random permutation
fil = fopen(sprintf('%s_%s.txt',mfilename,run_stamp),'a');
Dims = floor(n*(0.05:0.01:0.7));


fold = struct('train_ind',[],'T',[],'U',[],'K',[],'test_ind',[],'max_dim',[]);
pls_dim = max(Dims);
for i = 1:nfolds
    train_slct = ones(n,1);
    train_slct(floor((i-1)*n/nfolds)+1:floor(i*n/nfolds)) = 0;
    train_slct = logical(train_slct);

    fold(i).test_ind=rand_ind(~train_slct);
    fold(i).train_ind = rand_ind(train_slct);

    A_train_fold = A_train(fold(i).train_ind,:);
    train_lbl_fold = matrix_labels(fold(i).train_ind,:);

    disp 'time of KPLS'
    tic
    MAP = KPLS
    MAP = MAP.Learn(A_train_fold,train_lbl_fold,pls_dim);%,length(unique(letter_labels))-1);
    toc
    fold(i).MAP = MAP;
    %fold(i).U = MAP.U;
    fold(i).max_dim = size(MAP.T_,2);
    %fprintf(fil,'\n Norm(K)=%f',normK);
end

fold_dims = cell2mat({fold.max_dim});

min_dim = min(fold_dims);
Dims = Dims(Dims<=min_dim);

errors = zeros(1,length(Dims));

for j = round([2:9]\10*length(Dims))
    pls_dim = Dims(j);
    fprintf(fil,'\n\nFor pls_dim %d',pls_dim);
    for i = 1:nfolds
        disp 'time to reprepare matrices'
        tic
%         A_train_fold = A_train(fold(i).train_ind,:);
%         A_train_fold = bsxfun(@minus,A_train_fold,mean(A_train_fold));
% 
%         train_lbl_fold = matrix_labels(fold(i).train_ind,:);
%         train_lbl_fold = bsxfun(@minus,train_lbl_fold,mean(train_lbl_fold));

        A_test_fold = A_train(fold(i).test_ind,:);
        A_test_fold = bsxfun(@minus,A_test_fold,mean(A_test_fold));

        test_lbl_fold = labels_train(fold(i).test_ind);

        U_fold = fold(i).U(:,1:pls_dim);
        T_fold = fold(i).T(:,1:pls_dim);
        toc

        disp 'Reconstructing projection matrix M'
        tic
        M = A_train_fold'*(U_fold*inv(T_fold'*fold(i).K*U_fold))*(T_fold'*train_lbl_fold);
        toc
        fprintf(fil,'\n Norm(K)=%f',normK);

        %test using the test fold
        %test using classifying new_labels using max or using Ttrain and Ttest.

        est_test_lbl_fold= A_test_fold*M;
        [~,est_test_lbl_fold] = max(est_test_lbl_fold,[],2);
        diff = est_test_lbl_fold - test_lbl_fold;
        err_ = nnz(diff)/length(test_lbl_fold);

        errors(j) = errors(j)+err_;
    end
    errors(j)
    if j>20
        var(errors(j-7:j))
        %var(errors(1:j-7))
        if var(errors(j-7:j)) < v_thresh
            %errors(j-7:j)
            %min(errors)
            %if (errors(j)- min(errors))>m_thresh
                break;
            %end
        end
    end
end
% Choose the best_pls_dim
[~,best_dim_indx] = min(errors);
pls_dim = Dims(best_dim_indx);
fprintf(fil,'\n\n=========================== end of crossvalidation');
fprintf(fil,'\nBest pls_dim is %d',pls_dim);
fclose(fil);

    