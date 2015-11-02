function [A]=direct_gaussian_RBF(X,cent)
% Gaussian Radial Basis Function (Gaussian-RBF) interpolation
% The centres are assumed given in the n by dim array cent.
%
% Syntax  [A]=direct(cent)
%
% Arguments:         
%         X     m*dim array
%         cent  n*dim array  coordinates of centers
%               of reals     
% Output  A     m*n array for the radial basis function
%                            interpolant directly.
%
% Copyright (c) 2012 Amr M. Bakry

%r=sqrt(dist2(X,cent));
%A=exp(-0.5*(r.*r));
r=dist2(X,cent);
A=exp(-0.5*(r));
A= sparse(A);
