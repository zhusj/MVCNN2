function [newVec, dim] = shiftVec(X, dim)
% FUNCTION SHIFTVEC(X, DIM) This is to shift vector in unfold and fold for
% the shift of dimension vector
% by Chan-Su Lee

N = size(X,2);
newVec = [];
if dim>0
    newVec(1: N-dim) = X(dim+1:N);
    newVec(N-dim+1:N) = X(1:dim);
elseif dim<0
    newVec(1:-dim) = X(N+dim+1:N);
    newVec(-dim+1:N) = X(1:N+dim);
else
    newVec = X;
end