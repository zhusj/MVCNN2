function [foldedX] = fold(X, dim,sizeX)
% FUNCTION FOLD(X, DIM,SIZEX) fold X multidimensional array from 2Dimensional
% matrix unfolded in dim to sizeX multidimensional arrray
% - input:
%       + X: multidimensional array
%       + dim: base mode dimension used for row dimension
%       + sizeX: number of dimensions in multidimensional array
% - output:
%       + foldedX: folded new matrix
% by Chan-Su Lee
% - version 0.1: 02/09/04
%       + create primitive function based on shift function in Matlab
%


N = length(sizeX);
foldedX = [];

if N < dim
    disp('Error in folding. dimension to fold may not correct');
    return
end

% foldedX = reshape(shiftdim(reshape(X', shiftVec(sizeX,dim-1)),N-dim+1),sizeX);
foldedX = reshape(shiftdim(reshape(X, shiftVec(sizeX,dim-1)),N-dim+1),sizeX);