function storeFile = learn_KPLS(run_stamp,filename,dim,plot_flg,pls_dim)
% FILE NAME: HOFDA.m (Higher Order Fisher Diecreminating Analysis)
% tensor analysis for gcf not gscf
    diary 'run_log_learn_HOFDA.txt'
    echo on

    if ~exist('plot_flg','var')
        plot_flg = 1;
    end

    load(filename,'train_GCF');
    
    sz = size(train_GCF);
    %sz = sz(1:end-1);
    
    
    %prepare labels
    labels_train = get_unfolding_label(dim,sz);
    matrix_labels = get_matrix_label(labels_train,sz(dim));
    
    %prepare the data matrix
    B_train = List_CFs(train_GCF);
    
    clear train_GCF;
    %% Kernel computing
    n = size(B_train,1);
    %A_train=get_gram_kernel(B_train,B_train);%full(directs_grbf_reg(Y,Y));
    A_train=B_train;
    %A_train = bsxfun(@minus, A_train, mean(A_train));
    % r=dist2(X,X);
    % A = exp(-0.5*r);
    % clear r;
    %K_train = A_train*A_train';%get_gram_kernel(A,A);
    %K_train = cosine_similarity(A_train,A_train);%this is similar to %K_train = A_*A_' , where A_ = A_train after row norm normalizing
%     w = ones(n,1);
%     w= w*w';
%     w = eye(n)-w/n;
%     K_train = w*K_train*w';%normalizing K using equation (13) in the paper

    %% start KPLS processing
    [new_labels,MAP,K_train,normK] = KPLS(A_train,matrix_labels,pls_dim,1E-1);%,length(unique(letter_labels))-1);
    
    if plot_flg
        [~,lr] = max(new_labels,[],2);
        lr = mod(lr-1,26)+1;
        figure, plot(lr);
    end
    
%     letter_labels = mod(letter_labels-1,26)+1;
    
    storeFile = sprintf('Learn_KPLS_%s',run_stamp);
    save(storeFile,'MAP', 'B_train','A_train','K_train','labels_train','new_labels','normK');
    %% Plotting
    if plot_flg
        plotting(new_labels,labels_train);
    end
    echo off
    diary off
end
        
