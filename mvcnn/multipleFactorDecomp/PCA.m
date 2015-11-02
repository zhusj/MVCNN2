function final = PCA(X,OutputSize)
%[r, Columns] = size(x);  % find size of input matrix
m=mean(X);                  % find mean of input matrix
y=X-ones(size(X,1),1)*m;    % normalise by subtracting mean
c=cov(y);                   % find covariance matrix
%c = y*y';
[V,D]=eig(c);               % find eigenvectors (V) and eigenvalues (D) of covariance matrix
[D,idx] = sort(diag(D));    % sort eigenvalues in descending order by first diagonalising eigenvalue matrix, idx stores order to use when ordering eigenvectors
D = D(end:-1:1)';
V = V(:,idx(end:-1:1));     % put eigenvectors in order to correspond with eigenvalues
V2d=V(:,1:OutputSize);        % (significant Principal Components we use, OutputSize is input variable)
prefinal=V2d'*y';
final=prefinal';            % final is normalised data projected onto eigenspace
end