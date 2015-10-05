function x = normalizePoints(x)
% x - points (n points x m dimensions)
% normalize points (zero-mean and unit st. dev)
x = (x - repmat(mean(x),...
    [size(x,1) 1])) ./ repmat(std(x),[size(x,1) 1]);