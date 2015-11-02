
% learn GRBF mapping
% input:
% D : original data matrix Nxd  N: is the number of samples, d: is the dimensionality of the input images
% X : embedding space representation. Nxe matrix N: number of samples, e: the dimesnionality of the embedding space
% cent: centers: exm e: the dimesnionality of the embedding space

function CF=learnmapping_grbf_regulize(D,cent)
%CF = learnmapping_KPLS(D);
%return

%if centers number is not given then set it (better)
if nargin < 2
    ncent = 8;% mapping centers, Need to know if this number affects the performance???????
    icent=((1:ncent)-1)'*2*pi/ncent;
    cent = [cos(icent) sin(icent)];
    %Nb=length(cent);
    %d=2;
end

%%%%%
%cent = cent(2:end,:);
%%%%%
% % embed on a circle
% %measure how close every two consecutive frames:
d = dist2(D,D);
d = [diag(d,1);d(1,end)];
d = d/sum(d);
%accumlate distances so that last value is 1:
accumlate = tril(ones(length(d),length(d)),0);
d = accumlate*d;
d = [0;d(1:end-1)];

% frames=size(D,1);%angletrain{seq};
% d = (1:frames)'/frames;
%t=[1:frames]'*2*pi/frames;
%get embedding on circle.
t=d*2*pi;
X=[cos(t) sin(t)];


%implementation of higher dimensional regularization
%cf = (G'G+lambdag)-1 G'y
lambda = -50.0;
N=size(D,2);
G =directs_grbf_reg(X,cent);
g =directs_grbf_reg(cent,cent);
F =G'*D;
B =full(G'*G-lambda*g);
CF = B\F;
%%%%%
% CF = CF(2:end,:);
%CF = inv(G'*G-lambda*g)*F;
end
    
    

