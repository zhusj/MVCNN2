function S = cosine_similarity(X,Y)
%row wise data matrices X,Y. Data sample i is ith row.
nX = sqrt(sum(X.^2,2));
X = bsxfun(@rdivide, X,nX);
if exist('Y','var')
    nY = sqrt(sum(Y.^2,2));
    Y = bsxfun(@rdivide, Y,nY);
else
    Y = X;
end

S = (X*Y');