function [best_dim,result_list] = estimate_kpls_dim_cv(Data,Labels,split_labels,stp_pt,p_classification)
%function for  computing best PLS_dim using cross-validation

N= length(unique(split_labels));

kpls_list = cell(N,1);

min_err = inf;
max_dim = 100;
min_dim = 5;

err_list = zeros(max_dim ,1);

for i = 1:N-1
    validation_ind = (split_labels==i | split_labels==(i+1));
    train_ind = ~validation_ind;
    
    kpls_list{i} = KPLS;
    kpls_list{i} = kpls_list{i}.Learn(Data(train_ind,:),Labels(train_ind),max_dim,stp_pt,p_classification);

    % here i am comparing by how far we have the ability to predict the
    % labels
    new_max_dim = min(max_dim,kpls_list{i}.getDimension);
    new_min_dim = min(min_dim,new_max_dim);
    
    
    
    for dim = new_min_dim:new_max_dim
         kpls_list{i}.setDimension(dim);
        [~,estimate_valid_lbls] = kpls_list{i}.Predict(Data(validation_ind,:),'input');
        err = nnz(estimate_valid_lbls-Labels(validation_ind))/sum(validation_ind);
        if err<1E-12
            break;
        end
        err_list(dim) = err_list(dim)+ err;
        
        if err_list(dim)<min_err
            best_dim = dim;
            min_err = err_list(dim);
        end
    end
    
end
err_list(1:min_dim-1)=inf;
dim_regularize = exp(-(new_min_dim:new_max_dim)'/3);
err_list = err_list - dim_regularize;

[~,best_dim] = min(err_list);
result_list = err_list;