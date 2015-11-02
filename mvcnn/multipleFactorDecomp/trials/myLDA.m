function B = myLDA(A,L)
    uL = unique(L);
    W = zeros(length(uL),size(A,1),size(A,1));
    sz = zeros(length(uL));
    for i = 1:length(uL)
        r = A(:, L==uL(i));
        sz(i) = size(r,2);
        m = mean(r,2);
        m = repmat(m,1,sz(i));
        r = r - m;
        W(i,:,:) = (r*r');%/sz(i);
    end
    %sz = repmat(sz,[1,size(A,1),size(A,1)]);
    Sw = sum(W)/size(A,2);
    m = mean(A,2);
    m = repmat(m,1,size(A,2));
    r = A - m;
    Sb = (r*r')/size(A,2);
    
    
    [V,D] = eig(Sb,Sw);
    
    [lambda, ind] = sort(diag(D));
    lambda(isnan(lambda)) = 0;
    thresh = 1E-5;
%     indx = find(diag(lambda)>thresh);
% 	V = V(:,ind(min(indx)+1:end);%min([no_dims size(V, 2)])));    
    V = V(:,diag(lambda)>thresh);
    
	% Compute mapped data
	B = V'*A;
end
        
