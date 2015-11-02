function [newA] = tmul2(A, M, nmode)
% FUNCTION TMUL(A, M, NMODE) TENSOR MULTIPLIATION FOR GIVEN TENSOR A and
% MATRIX M WITH RESPECT NMODE ELEMENT
% by Chan-Su Lee

sizeA = size(A);

length = size(M,2);

newA = [];
tempA = [];
if length < nmode
    disp('Error in tmul2: The dimension of mode array and mode should be same')
    return
end

tempA = A;
for i = 1:length
    if i ~= nmode
        tempA = tmul(tempA,M{i},i);
    end
end

[newA] = svdn(tempA,nmode,size(M{nmode},2));

