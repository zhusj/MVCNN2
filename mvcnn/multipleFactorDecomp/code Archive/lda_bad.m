function [mappedX, mapping] = lda(X, labels, no_dims)
%We follow the method mensioned in paper A Direct LDA algorithm for
%high-deimensional data - with application to face recognition
   
    %handle_traspose = 0;
% 	if size(X,1)<size(X,2)
%         handle_traspose = 1;
%         X_orig = X;
%         X = X';
%     end
   
    
	% Make sure data is zero mean
    mapping.mean = mean(X, 1);
	X = bsxfun(@minus, X, mapping.mean);
	X_old = X;
    X = directs_grbf_reg(X,X);
	% Make sure labels are nice
	[classes, ~, labels] = unique(labels);
    nc = length(classes);
    
    if ~exist('no_dims', 'var') || isempty(no_dims)
        no_dims = nc-1; % output space dimension should be < nc
    end
	% Intialize Sw
	%Sw = zeros(size(X, 2), size(X, 2));
    
    % Compute total covariance matrix
    %St = cov(X);
    means = zeros(nc,size(X,2));
    
    X2 = zeros(size(X));
	% Sum over classes
    for i=1:nc
        % Get all instances with class i
        cur_X = X(labels == i,:);
        means(i,:) = sqrt(size(cur_X,1))*mean(cur_X);

        X2(labels == i,:) = cur_X - repmat(means(i,:),size(cur_X,1),1);
        means(i,:) = means(i,:);% - mapping.mean;
        % Update within-class scatter
        %C = cov(cur_X);
        %p = size(cur_X, 1) / (length(labels) - 1);
        %Sw = Sw + (p * C);
    end
    %Sw is m x m
    
    %as the feature dimensionality is very hight, so we compute the eig of
    %Sb as follows:
    %basically, Sb = means'*means
    %so we can compute eig of Sb as
    Sb_ = means*means';
    [V_,D_] = eig(Sb_); %V_ is nc X nc
%     [~, ind] = sort(diag(D_), 'descend');
%     %V_ = V_';
% 	V_ = V_(:,ind(1:no_dims));
%     D_ = diag(D_);
%     Db = diag(D_(1:no_dims));
%     diag(Db);
   Db = D_; 
    Vb = V_'*means; %V is nc x m
        
    Z = (Db^(-0.5))*Vb; % Z is nc X m
    
    %test Z*Sb*Z' = I
    %Z * (means'*means) * Z'
    
    %Sw_  = Z * Sw*Z'; %Sw_ is nc x nc
    X3 = Z*X2'; % X3 is nc x n
    Sw_  = X3*X3'; %Sw_ is nc x nc
    [V_,D_] = eig(Sw_); % V_ is nc x nc
    %V_ = V_';
    %exclude one
    D_(isnan(D_)) = 0;
	[~, ind] = sort(diag(D_), 'descend');
	V_ = V_(:,ind(1:no_dims));
    
    A = V_'* Z; % A is nc x m
    
    
    mapping.M = real(A);
    mappedX = real(A)*X_old;
    
end