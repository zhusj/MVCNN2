function [CoreMatrix, modebasis,modelsv] = HOSVD(D, sizeD)
% FUNCTION HOSVD(D) HIGHER-ORDER SINGULAR VALUE DECOMPOSITION FOR GIVEN
% MULTILINEAR MATRIX D
% - input
%      + D:multilinear matrix
% - output
%      + CoreMatrix: CoreMatrix in HOSVD
%      + modebasis: cell array has basis for each mode
% - version 0.1: 02/11/04
%      + initial implementation based on Vasilescu paper
% by Chan-Su Lee


if nargin <2
    sizeD = size(D);
end

modebasis = {};
modelsv = {};

Dlength = length(sizeD);

for i=1:Dlength
    [modebasis{i}, modelsv{i}] = svdn(D,i);
end

CoreMatrix = D;
for i=1:Dlength
    CoreMatrix = tmul(CoreMatrix,modebasis{i},i);
end

    