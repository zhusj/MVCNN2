% learn Kernel PLS mapping: implementatin to what in paper Kenrk Partial
% Least square regression in RKHS.
% input:
% D : original data matrix Nxd  N: is the number of samples, d: is the dimensionality of the input images
% X : embedding space representation. Nxe matrix N: number of samples, e: the dimesnionality of the embedding space
% cent: centers: exm e: the dimesnionality of the embedding space

function CF=learnmapping_KPLS(Y)
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
d = dist2(Y,Y);
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


%% normalize
X = bsxfun(@minus, X, mean(X));
Y = bsxfun(@minus, Y, mean(Y));
%% estimate 
n = size(X,1);
A=full(directs_grbf_reg(X,cent));
thresh = 1E-15;
t_old = zeros(n,1);
t = ones(n,1);
u = rand(n,1);
T = [];
U = u;
u_old = -1*u;
indx = 1;
K = A*A';
w = ones(n,1);
w= w*w';
w = eye(n)-w/n;
K = w*K*w';
while norm(t-t_old)>thresh || norm(u-u_old)>thresh
    t_old = t;
	t = K*u;
	t = t/norm(t);
    T = [T,t];
	c = Y'*t;
    u_old = u;
	u = Y*c;
	u = u/norm(u);
    U = [U,u];
    indx = indx +1;
end
U = U(:,1:end-1);
v = (eye(n)-t*t');
K = v*K*v;
Y = v*Y;

B = A'*U*inv(T'*K*U)*T'*Y;
CF = B;

norm(Y-A*B)
end