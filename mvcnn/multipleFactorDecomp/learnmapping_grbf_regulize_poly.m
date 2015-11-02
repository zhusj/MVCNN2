
% learn GRBF mapping
% input:
% D : original data matrix Nxd  N: is the number of samples, d: is the dimensionality of the input images
% X : embedding space representation. Nxe matrix N: number of samples, e: the dimesnionality of the embedding space
% cent: centers: exm e: the dimesnionality of the embedding space

function [CF,new_data]=learnmapping_grbf_regulize_poly(Data,cent,points,lambda)
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
    ncent = 7;% mapping centers, Need to know if this number affects the performance???????
    icent=((1:ncent)-1)'*2*pi/ncent;
    cent = [cos(icent) sin(icent)];
    %Nb=length(cent);
    %d=2;
end

% Set embedding on a unit circle
if ~exist('points','var')||isempty(points)
    frames=size(Data,1);
    d = ((1:frames)'-1)/frames;
    %get embedding on circle.
    t=d*2*pi;
else
    t = points;
end

% %make sure that the input points temporally sorted
% [t,ord] = sort(t);
% Data = Data(ord,:);
% X=[cos(t) sin(t)];
X = t;

%implementation of higher dimensional regularization
% paper: T. Poggio and F. Girosi. Networks for approximation and 
% learning. Proceedings of the IEEE, pages 1481–1497, 1990.
%cf = inv(G'G+lambdag) G'y
%initialize the regularizer
%lambda = 20;

N=size(Data,2);
dim = size(cent,2);
%learn nonlinear mapping
%for different vlues of regularizer lambda
A=directs_grbf_reg_poly(X,cent,lambda);

F=[Data; zeros(dim+1,N)];
%CF*A = F
Af=full(A);
CF=Af\F;

new_data = Af*CF;
new_data = new_data(1:end-3,:);
end
    
    

