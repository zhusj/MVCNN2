%direct.m - Forming the interpolation matrix

function [A]=directs_grbf_jacobian(X,cent)

% get the jacobian of the rbf kernel with basis cents at point X

% Form the (n+dim+1)*(n+dim+1) matrix corresponding to 
% RBF interpolation with basic function phi and linears.
%
% The centres are assumed given in the n by dim array cent.
% phi is assumed given as a function of r. It is coded in 
% the Matlab function phi. m
%
% Syntax  [A]=direct(cent)
%
% Input          
%         cent  n*dim array  coordinates of centers
%               of reals     
% Output  A     (n+dim+1)*   Symmetric matrix of the linear
%               (n+dim+1)    system obtained if we solve
%               array of     for the radial basis function
%                            interpolant directly.
%
% Write the matrix A in the form
%             B    P
%     A   =  
%             P^t  O
% where P is the polynomial bit.
%
[n dim]=size(cent);
m=size(X,1);

A=zeros(m,n);

r=sqrt(dist2(X,cent));
A=phi(r);

Xv = X'; Xv = Xv(:);
C = repmat(cent',m,1);

diff = bsxfun(@minus,C,Xv); %sigma^2 = 2, otherwise constant has to be verified
%a = A';a = a(:);
A = imresize(A,[m*dim,n],'nearest');%%@@TODO: REVALIDATE IN CASE OF MULTIPLE POINTS
A = A.*diff;
%A = bsxfun(@times,diff,a);
%A= sparse(A);

%A = A + lambda*eye(size(A));

% %
% % Now the polynomial part
% %
%  P1=[ones(m,1) X];
%  P2=[ones(n,1) cent];
%  A = [ sparse(A)      sparse(P1);
%        sparse(P2')  sparse(dim+1,dim+1)];


