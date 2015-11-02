
% learn GRBF mapping
% input:
% D : original data matrix Nxd  N: is the number of samples, d: is the dimensionality of the input images
% X : embedding space representation. Nxe matrix N: number of samples, e: the dimesnionality of the embedding space
% cent: centers: exm e: the dimesnionality of the embedding space

function CF=learnmapping_grbf_regulize_poly(D,cent)
%if centers number is not given then set it (better)
if nargin < 2
    ncent = 10;% mapping centers, Need to know if this number affects the performance???????
    icent=((1:ncent)-1)'*2*pi/ncent;
    cent = [cos(icent) sin(icent)];
    %Nb=length(cent);
    %d=2;
end


% embed on a circle
%measure how close every two consecutive frames:
d = dist2(D,D);
d = [diag(d,1);d(1,end)];
d = d/sum(d);
%accumlate distances so that last value is 1:
accumlate = tril(ones(length(d),length(d)),0);
d = accumlate*d;
d = [0;d(1:end-1)];
%frames=size(D,2);%angletrain{seq};
%t=[1:frames]'*2*pi/frames;
%get embedding on circle.
t=d*2*pi;
X=[cos(t) sin(t)];

%learn nonlinear mapping
%for different vlues of regularizer lambda
error_min= inf;
opt_lambda = [];
opt_CF = [];
for lambda = -2.0%0:0.1:1
    [n dim] =size(cent);
    N=size(D,2);
    A=directs_grbf_reg_poly(X,cent,lambda);

    F=[D; zeros(dim+1,N)];
    %CF*A = F
    Af=full(A);
    CF=Af\F;
    rec_err = norm(F-Af*CF)
    if rec_err < error_min
        error_min = rec_err;
        opt_lambda = lambda;
        opt_CF = CF;
    end
end
%opt_lambda
%error_min
CF = opt_CF;
end
    
    

