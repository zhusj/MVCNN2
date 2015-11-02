function [CoreMatrix, modebasis,modelsv] = HOSVDDR2(D, sizeD, count_num)
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

if nargin < 3
    count_num = 5;
end

count = 1;

if nargin <2
    sizeD = size(D);
end

Dlength = length(sizeD);
modebasis = cell(1,Dlength);
modelsv = cell(1,Dlength);

for i=1:Dlength
    [modebasis{i}, modelsv{i}] = svdn(D,i,sizeD(i));
end


flag = ones(1,Dlength);

imodebasis = modebasis;
%imodelsv = modelsv;
if nargin  >= 2  % in case of dimensionality reduction
    while count <= count_num
        for i=1:Dlength
            %sprintf('HOSVDDR %d iteration',i)
            imodebasis{i} = tmul2(D, modebasis, i);
            %diff = sum(abs(diag(imodebasis{i}'*modebasis{i})))
        end % for  
        modebasis = imodebasis;
        %for i=1:Dlength
        %    modebasis{i} = imodebasis{i};    
        %end
        count = count +1
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

d = CoreMatrix;
for i=1:Dlength
    d = tmul(d,modebasis{i}',i);
end
disp 'Tensor reconstruction error is:'
norm(unfold(d,1)-unfold(D,1))
    