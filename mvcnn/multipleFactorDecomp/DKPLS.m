% learn Kernel PLS mapping: implementatin to what in paper Kenrk Partial
% Least square regression in RKHS.
% input:
% D : original data matrix Nxd  N: is the number of samples, d: is the dimensionality of the input images
% X : embedding space representation. Nxe matrix N: number of samples, e: the dimesnionality of the embedding space
% cent: centers: exm e: the dimesnionality of the embedding space

function [Y,B]=DKPLS(A,K,L,D)
%% normalize
B.mean = mean(L);
L = bsxfun(@minus,L, B.mean);
B.sgma = sqrt(var(L));
%L = L/B.sgma;
L =  bsxfun(@rdivide,L,B.sgma);
%X = bsxfun(@minus, X, mean(X));
%Y = bsxfun(@minus, Y, mean(Y));

% w = r'*y;
% w = w / (w'*w);
%r_new = r*(w*w');
% W = X'*L;
% W = W/(W'*W);
% B.M = W*W';
% Y = X*B.M;
% return;
if ~exist('D','var') || (D>min(size(A)))
    D = min(size(A));
end
n = size(A,1);
%% Iterating
K1 = K;
L1 = L;
thresh = 1E-10;
indx = 1;
T = [];
U = [];
flag=0;
while indx<=D
%     t_old = zeros(n,1);
%     t = ones(n,1);
%     u = rand(n,1);
%     u_old = -1*u;

%     while norm(t-t_old)>thresh || norm(u-u_old)>thresh
%         t_old = t;
        t = K*K'*L;
        if norm(t)<thresh
            flag = 1;
            break;
        end
        t = t/norm(t);
        v = (eye(n)-t*t');
        K = v*K;
        %c = L'*t;
        u_old = u;
        u = v*L;
        if norm(u)<thresh
            flag = 1;
            break;
        end
        u = u/norm(u);
    end
    if flag
        D = indx-1
        break;
    end
    U = [U,u];
    T = [T,t];
    indx = indx +1;
    
    L = v*L;
end
%U = U(:,1:end-1);
B.M2 = U*inv(T'*K1*U)*T'*L1;
B.M = A'*B.M2;
Y = A*B.M;
%B.mean = zeros(1,size(Y,2));

%norm(Y-A*B);
end