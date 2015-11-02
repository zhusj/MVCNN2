function [unfoldedX, sizeX, dim] = unfold(X, dim)
% FUNCTION FOLD(X, DIM) fold X multidimensional array to 2Dimensional
% matrix based on given dim
% - input:
%       + X: multidimensional array
%       + dim: base mode dimension used for row dimension
% - output:
%       + unfoldedX: unfolded new matrix
%       + sizeX: original multidimensional size, which can be used for
%       folding
%       + dim: base mode dimension
% by Chan-Su Lee
% - version 0.1: 02/09/04
%       + create primitive function based on shift function in Matlab
%


sizeX = size(X);

N = length(sizeX);
unfoldedX = [];

if N < dim
    disp('Error in unfolding. dimension to unfold may not correct');
    return
end

% unfoldedX = reshape(shiftdim(X, dim-1), prod(sizeX)/sizeX(dim),sizeX(dim))';
% sizeX = shiftVec(sizeX,dim-1);

unfoldedX = reshape(shiftdim(X, dim-1), sizeX(dim), prod(sizeX)/sizeX(dim));
sizeX = shiftVec(sizeX,dim-1);