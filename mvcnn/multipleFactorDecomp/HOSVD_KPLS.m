function [CoreMatrix, modebasis,modelsv] = HOSVD_KPLS(D, sizeD, modeLabel,count_num)
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

if ~exist('count_num','var')
    count_num = 5;
end

count = 1;
dimenReduct = 1;
if ~exist('sizeD','var')
    sizeD = size(D);
    dimenReduct = 0;
end
%dimenReduct = 0;%%%%%should be removed
Dlength = length(sizeD);
modebasis = cell(1,Dlength);
modelsv = cell(1,Dlength);

DS = D;

%while count < count_num
for i=1:Dlength
    if ~isempty(modeLabel{i})
        [modebasis{i}] = KPLSn(DS,i, modeLabel{i}, sizeD(i));
        DS = tmul(DS,modebasis{i},i);
    else
        [modebasis{i}] = svdn(DS,i,sizeD(i));
    end
end
%count = count +1;
%end

flag = ones(1,Dlength);

imodebasis = modebasis;
%imodelsv = modelsv;
if dimenReduct  % in case of dimensionality reduction
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
        %plotting(modebasis{1},[1:26,1:26]',0);
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