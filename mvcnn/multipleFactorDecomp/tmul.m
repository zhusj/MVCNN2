function [newA, newSizeA] = tmul(A, M, nmode)
% FUNCTION TMUL(A, M, NMODE) TENSOR MULTIPLIATION FOR GIVEN TENSOR A and
% MATRIX M WITH RESPECT NMODE ELEMENT
% by Chan-Su Lee

sizeA = size(A);

[mj, mi] = size(M);
newA = [];
% if sizeA(nmode) ~= mi
%     disp('Error in tmul: The dimension of matrix and tensor in the mode does not same')
%     return
% end
% 

% mulA = M*unfold(A,nmode);
% changed for tmul2-verify which one is correct
mulA = M'*unfold(A,nmode);
newSizeA = sizeA;
newSizeA(nmode) = mi;

newA = fold(mulA, nmode, newSizeA);

