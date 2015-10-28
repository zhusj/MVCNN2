
function [X]=CtoX(CF,npoints,cent)
% Compute the drivative of the Gaussian RBF parameterization with respect
% to the input images
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

G =direct_gaussian_RBF(t,cent);
X = G*CF;

end
