function A = gaussian_rbf(X,cent)

d = size(cent,1)
n = size(X,1)
if size(X,2)>1
    r = sqrt(pdist2(X,cent,'cosine'));
else
    r = abs(bsxfun(@minus,X,cent'));
end

if d==n
    W = eye(d);
else
    thresh = 1E-15;
    
    Q = r*r';
    [U,S,V] = svd(Q);
    s = diag(S);
    l = min(find(s<thresh))
    if isempty(l)
        W = V(:,end-d+1:end);
    else
        W = V(:,l-d:l);
    end
end

M = W*W';
A = sparse(exp(r'*M*r));



