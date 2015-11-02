% learn Kernel PLS mapping: implementatin to what in paper Kenrk Partial
% Least square regression in RKHS.
% input:
% D : original data matrix Nxd  N: is the number of samples, d: is the dimensionality of the input images
% X : embedding space representation. Nxe matrix N: number of samples, e: the dimesnionality of the embedding space
% cent: centers: exm e: the dimesnionality of the embedding space

function [Y,MAP,K1,normK]=KPLS(A,L,D,stp)
%% normalize
%X = bsxfun(@minus, X, mean(X));
%Y = bsxfun(@minus, Y, mean(Y));
%B = [];
% w = r'*y;
% w = w / (w'*w);
%r_new = r*(w*w');
% W = X'*L;
% W = W/(W'*W);
% B.M = W*W';
% Y = X*B.M;
% return;

if ~exist('D','var') || isempty(D) || D<1 || (D>=min(size(A))) 
    D = min(size(A));
else
    D = D+1;
end

if ~exist('stp','var') || stp==0 || isempty(stp)
    stp = 1E-4;
end
%D = 101;
n = size(A,1);

%data normalization
mn = mean(L);
L = bsxfun(@minus,L, mn);
%sgma = sqrt(var(L));
%L = L/B.sgma;
%L =  bsxfun(@rdivide,L,sgma);

%K = A*A';
K = get_gram_kernel(A);
w = ones(n,1);
w= w*w';
w = eye(n)-w/n;
K = w*K*w';%normalizing K using equation (13) in the paper
%% Iterating
%A1 = A;
A1 = bsxfun(@minus, A, mean(A));
K1 = K;
L1 = L;
thresh = 1E-8;
indx = 1;
T = [];
U = [];
C = [];
%W = [];
S = [];
flag=0;
reason = 1; % indicate the reason of stopping the iteration.
while indx<D
    t_old = zeros(n,1);
    t = ones(n,1);
    u = rand(n,1);
    u_old = -1*u;
    while norm(abs(t)-abs(t_old))>thresh %|| norm(u-u_old)>thresh
        t_old = t;
        %w = A'*u;
        t = K*u;
        if norm(t)<thresh
            flag = 1;
            reason = 2;
            break;
        end
        s = norm(t);
        t = t/s;
        
        c = L'*t;
        u_old = u;
        u = L*c;
%         if norm(u)<thresh
%             flag = 1;
%             break;
%         end
        u = u/norm(u);
    end
    if flag
        D = indx-1
        break;
    end
    
    U = [U,u];
    T = [T,t];
    S = [S,s];
    C = [C,c];
    %w = A'*u_old;
    %w = w/norm(w);
    %W = [W,w];
    
   
    v = (eye(n)-t*t');
    K = v*K*v;
    %A = v*A;
    L = v*L;
    normK = norm(K);
    if normK<stp%thresh
%         %flag = 1;
        reason = 3;
        D = indx;
        break;
    end
    indx = indx +1;
end
%%% get the last dimension - not good this diemnsion is not completely orthogonal to
%%% others
% u = rand(n,1);
% t = K*u;
% t = t/norm(t);
% c = L'*t;
% u = L*c;
% u = u/norm(u);
% t = K*u;
% t = t/norm(t);
% U = [U,u];
% T = [T,t];
% C = [C,c];

disp 'Reason of exit loop is'
reason

disp 'Number of iteration:'
D
%U = U(:,1:end-1);
MAP.R2 = U*inv(T'*K1*U);
MAP.M2 = U*inv(T'*K1*U)*T'*L1;
if numel(size(A1))==2 %if A1 is not tensor
    MAP.R = A1'*MAP.R2;
    MAP.M = A1'*MAP.M2;
else
    MAP.R = [];
    MAP.M = [];
end
MAP.T = T;
MAP.U = U;

%Y = A*B.M;
%R = A'*U*inv(T'*K1*U);
Y = K1*MAP.M2;
%W = A'*T;
%B.mean = zeros(1,size(Y,2));

%norm(Y-A*B);
end