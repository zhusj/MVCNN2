function A = vectorize_CF(GCF)
    sz = size(GCF);
    data_sz = sz(end);
    A = reshape(GCF,[],data_sz);
    %A_train=get_gram_kernel(B_train,B_train);%full(directs_grbf_reg(Y,Y));
    
    %A_train = bsxfun(@minus, A_train, mean(A_train));
    % r=dist2(X,X);
    % A = exp(-0.5*r);
    % clear r;