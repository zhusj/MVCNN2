function [U,S,V] = svdn(D, mode,dimension)
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

THRESHOLD_DIM = 200;

sizeD = size(D);
lengthD = length(sizeD);

U = [];
S = [];
V = [];
if nargin <3
    dimension = 0;
end

if lengthD < mode
    disp('Error in svdn: in correct mode index');
    return;
end

unfoldedD = unfold(D, mode);


% [U, S, V] = svd(unfoldedD,0);
%%%%% implement based the the size of ni, nj to be economic %%%%%%%

% if dimension ~= 0 % optimization for lower dimension
if dimension ~= 0 % optimization for lower dimension
    [ni, nj] = size(unfoldedD);
    if ni >nj
        DTD = unfoldedD'*unfoldedD;
%         [V, SS, VT] = svd(DTD);
        [V, SS, VT] = svds(DTD,dimension);
        S = sqrt(SS);
        invS = inv(S);
        U = unfoldedD*V*invS;
    else 
        DDT = unfoldedD*unfoldedD';
%         [U, SS, UT] = svd(DDT);
        [U, SS, UT] = svds(DDT, dimension);
        S = sqrt(SS);
        invS = inv(S);
        V = invS*U'*unfoldedD;
    end
else
%     DDT = unfoldedD*unfoldedD';
%     [U, SS, UT] = svd(DDT,0);

    [U, S, V] = svds(unfoldedD,dimension);
%     [U, S, V] = svd(unfoldedD,dimension);


%     S = sqrt(SS);
%     invS = inv(S);
%     V = invS*U'*unfoldedD;
end
