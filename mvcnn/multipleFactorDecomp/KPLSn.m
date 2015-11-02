function [T,U,W,B] = KPLSn(D, mode, modeLabel, dimension)
% FUNCTION SVDN N-MODE SINGULAR VALUE DECOMPOSITION
% - input
%      + D: multilinear matrix
%      + mode: dimension to unfold for SVD
% - output
%      + U: basis for the given mode
%      + S: diagonal term
%      + V: V in SVD
% - version 0.1: 02/11/04
%      + initial implementation for HOSVD
% - version 0.2: 02/23/04
%      + add option with reduced dimensionality
% by Chan-Su Lee

sizeD = size(D);
lengthD = length(sizeD);

U = [];
W = [];
T = [];
B = [];

if lengthD < mode
    disp('Error in svdn: in correct mode index');
    return;
end

A = unfold(D, mode);
[r c] = size(A);

%if r<=c
    n = r;%min(size(A));
    %K = A*A';%get_gram_kernel(A,A);
% else
%     DTD = A'*A;
% %         [V, SS, VT] = svd(DTD);
%     [V, SS, VT] = svds(DTD,dimension);
%     S = sqrt(SS);
%     invS = inv(S);
%     T = A*V*invS;
%     return;
% end

if ~exist('dimension','var') || dimension>n
    dimension = n;
end


% w = ones(n,1);
% w= w*w';
% w = eye(n)-w/n;
% K = w*K*w';%normalizing K using equation (13) in the paper
[lbl_new,B,T,U] = KPLS(A,modeLabel, dimension);


