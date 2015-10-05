function [err] = manifold_alignment(x, theta, metric) 
% example call: [error] = manifold_alignment(pts, theta, 'correlation') 
% x - points (n points x m dimensions)
% angles - n angles x 1 - in radians
% metric - 'correlation' (only this is supported now)

% normalize points (zero-mean and unit st. dev)
x = (x - repmat(mean(x),...
    [size(x,1) 1])) ./ repmat(std(x),[size(x,1) 1]);

% sample pts on a unit-circle
theta = linspace(0,2*pi,12)';
pts_unit_circle = [cos(theta), sin(theta)];
pts_unit_circle = (pts_unit_circle - repmat(mean(pts_unit_circle),...
    [size(pts_unit_circle,1) 1])) ./ repmat(std(pts_unit_circle),[size(pts_unit_circle,1) 1]);
kernel_pts_unit_circle = pdist2(pts_unit_circle,pts_unit_circle);
kernel_pts_unit_circle = exp(-kernel_pts_unit_circle ./ (2*mean(kernel_pts_unit_circle(:))));

switch metric
    case 'correlation'
%         kernel_x = sqdistance(x',x');
        kernel_x = pdist2(x,x);
        kernel_x = exp(-kernel_x ./ (2*mean(kernel_x(:))));
        alignment = corr2(kernel_pts_unit_circle,kernel_x);
        err = 1 - alignment;    
    otherwise
        disp('this metric does not exist');
end