
function [dCF]=dCdLambda(CF,lambda,npoints,cent)
% Learn Gaussian RBF mapping with regularization
% Paper: 
%    A. Bakry; A. Elgammal, 
%    "MKPLS: Manifold Kernel Partial Least Squares for Lipreading and Speaker Identification,"
%    Computer Vision and Pattern Recognition (CVPR), 2013 IEEE Conference on , 
%    vol., no., pp.684,691, 23-28 June 2013

% Arguments:
%  Data : original data matrix Nxd  
%  points : embedding space representation. Nxe matrix
%  cent: centers: exm e: the dimesnionality of the embedding space
% Where
%  N: is the number of points (frames in our case), 
%  d: is the dimensionality of the image feature space
%  e: the dimesnionality of the embedding space
%
% Copyright (c) 2012 Amr M. Bakry

%if centers number is not given then set it (better)
if ~exist('cent','var') || isempty(cent)
    ncent = 12;% mapping centers, Need to know if this number affects the performance???????
    icent=((1:ncent)-1)'*2*pi/ncent;
    cent = [cos(icent) sin(icent)];
    %Nb=length(cent);
    %d=2;
end

% Set embedding on a unit circle
if ~exist('npoints','var')||isempty(npoints)
    npoints = 12;
end
d = ((1:npoints)'-1)/npoints;
%get embedding on circle.
t=d*2*pi;
t = [cos(t) sin(t)];

%implementation of higher dimensional regularization
% paper: T. Poggio and F. Girosi. Networks for approximation and 
% learning. Proceedings of the IEEE, pages 1481–1497, 1990.
%cf = inv(G'G+lambdag) G'y
%initialize the regularizer
if ~exist('lambda','var')
    lambda = .5;
end

G =direct_gaussian_RBF(t,cent);
g =direct_gaussian_RBF(cent,cent);

%F =G'*X;
B =full(G'*G+lambda*g);
B = B + 1E-5*eye(size(g));%add tiny diagonal to increase the stability of the system
%CF = B\F;
P = B\g;
dCF = P*CF;
%new_data = G*CF;

end
