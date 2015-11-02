function [CoreMatrix, modebasis,modelsv] = HOSVDDR(D, sizeD)
% FUNCTION HOSVD(D) HIGHER-ORDER SINGULAR VALUE DECOMPOSITION WITH DIMENSIONALITY REDUCTION FOR GIVEN
% MULTILINEAR MATRIX D
% - input
%      + D:multilinear matrix
% - output
%      + CoreMatrix: CoreMatrix in HOSVD
%      + modebasis: cell array has basis for each mode
% - version 0.1: 02/11/04
%      + initial implementation based on Vasilescu paper
% - version 0.2: 02/23/04
%      + find reduced dimension given by sizeD dimension
% - version 0.3: 08/18/04
%      + changed from svd to svdn to count reduced dimension in SVD
% by Chan-Su Lee

EPSILON = 0.002;
MAX_COUNT = 50;

count = 1;
if nargin <2
    sizeD = size(D);
end

modebasis = {};
modelsv = {};

Dlength = length(sizeD);

for i=1:Dlength
%     i
%     sizeD(i)
    [modebasis{i}, modelsv{i}] = svdn(D,i,sizeD(i));
end


flag = ones(1,Dlength);

imodebasis = {};
imodelsv = {};
if nargin  >= 2  % in case of dimensionality reduction
    while sum(flag) > 0 && count < MAX_COUNT
        for i=1:Dlength
            sprintf('HOSVDDR %d iteration',i)
            if flag(1,i) > 0
%                 [imodebasis{i}] = tmul2(D, modebasis, i)
                imodebasis{i} = tmul2(D, modebasis, i)
                %%% reconstruct the original higher order tensor and find
                %%% difference in the norm to find out accuracy compared to
                %%% the original ones
                diff = sum(abs(diag(imodebasis{i}'*modebasis{i})));
                if diff > (1-EPSILON)*sizeD(i)
                    flag(1,i) = 0;
                end
            end % if flag
        end % for  
        
        for i=1:Dlength
            modebasis{i} = imodebasis{i};    
        end
        count = count +1;
    end % while
% else
% 	for i=1:Dlength
%         [modebasis{i}, modelsv{i}] = svdn(D,i,sizeD(i));
% 	end    
end

CoreMatrix = D;
for i=1:Dlength
    CoreMatrix = tmul(CoreMatrix,modebasis{i},i);
end

    