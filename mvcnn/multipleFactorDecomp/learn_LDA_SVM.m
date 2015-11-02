function storeFile = learn_LDA_SVM(run_stamp,filename,learn_dim,pls_dim,p_cross_validation,split_dim_nm,plot_flg)
% FILE NAME: HOFDA.m (Higher Order Fisher Diecreminating Analysis)
% tensor analysis for gcf not gscf
global database
    diary 'run_log_learn_HOFDA.txt'
    echo on

    if ~exist('p_cross_validation','var')
        p_cross_validation = 0;
    end
    
    if ~exist('plot_flg','var')
        plot_flg = 0;
    end
    train_SVM = 0;
    
    
    storeFile = sprintf('%s_%s',mfilename, run_stamp);
    
    load(filename,'train_GCF');
    
    sz = size(train_GCF);
    sz = sz(1:end-1);
    
    
    %prepare labels
    if numel(size(train_GCF)) > 2
        train_labels = get_unfolding_label(learn_dim,sz);
    else
        load(filename,'train_labels');
    end
    
    matrix_labels = get_matrix_label(train_labels);
    
    if p_cross_validation && exist('split_dim_nm','var')
        split_dim = strmatch(split_dim_nm,{database.dimensions.name});
        split_labels = get_unfolding_label(split_dim,sz);
    end
    %prepare the data matrix
    A_train = list_CFs(train_GCF);
    
    clear train_GCF;
    %% Kernel computing
    n = size(A_train,1);
    
    %K_train = A_train*A_train';%get_gram_kernel(A,A);
    %K_train = cosine_similarity(A_train,A_train);%this is similar to %K_train = A_*A_' , where A_ = A_train after row norm normalizing
%     w = ones(n,1);
%     w= w*w';
%     w = eye(n)-w/n;
%     K_train = w*K_train*w';%normalizing K using equation (13) in the paper

    %% start KPLS processing
    if p_cross_validation
        pls_dim = get_PLS_dim_CV(A_train,train_labels,split_labels);
    end
    %[new_labels,MAP] = lda(A_train,matrix_labels,pls_dim);%,length(unique(letter_labels))-1);
    [new_labels,MAP] = lda(A_train,train_labels,pls_dim);%,length(unique(letter_labels))-1);
    MAP.T = new_labels;
    %fprintf(fil,'\n Norm(K)=%f',normK);
    
    pls_dim = size(MAP.T,2);
    
    save(storeFile,'MAP', 'train_labels','new_labels','pls_dim');
    
    if train_SVM
        %choose the best C_SVM parameters
        bestc=1; bestg=10;
        fprintf('\n\nSVM grid parameter search:');
        pls_dim = size(MAP.T,2);
        est_g_value = mean(var(MAP.T))*pls_dim; %mean variance of features*number of features;
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


        %SVM_model = ovrtrain(matrix_labels, MAP.T,sprintf('-c %g -g %g',bestc, bestg));
        SVM_model1 = ovrtrain(train_labels, MAP.T,sprintf('-c %g -g %g -q',bestc, bestg));
        %SVM_model2 = ovrtrain(labels_train, MAP.T,sprintf('-g 18'));

        if plot_flg
            [~,lr] = max(new_labels,[],2);
            lr = mod(lr-1,26)+1;
            figure, plot(lr);
        end

    %     letter_labels = mod(letter_labels-1,26)+1;



        save(storeFile,'SVM_model1','bestc', 'bestg','-append');
    end
    %% Plotting
    if plot_flg
        plotting(new_labels,train_labels);
    end
    echo off
    diary off
end
        
