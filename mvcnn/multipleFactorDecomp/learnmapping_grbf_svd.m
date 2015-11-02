
function CF=learnmapping_grbf_svd(D,X,cent)
%if centers number is not given then set it (better)
if nargin < 3
    ncent = 10;% mapping centers, Need to know if this number affects the performance???????
    icent=[1:ncent]'*2*pi/ncent;
    cent = [cos(icent) sin(icent)];
    %Nb=length(cent);
    %d=2;
end
% learn nonlinear mapping

[n dim] =size(cent);

N=size(D,2);

CF=[];

A=directs_grbf(X,cent);

%F=[D; zeros(dim+1,N)];

F=[D; zeros(dim+1,N)];

Af=full(A);

[U,S,V] = svd(Af,0);

s1=diag(S);

Nr=max(find(s1>1e-4));


S1=sparse([1:Nr]',[1:Nr]',(1./s1(1:Nr)),n+dim+1,n+dim+1);



CF=(U*S1*V')'*F;