function data_fit_svd
    diary 'results-SVD.txt'
    disp('***Fitting data using SVD method.');
    %Set up data point
    format longeng
    m = 51;
    x = .02*([1:m]'-1);
    y = 1+x;
    %plot(x,y)
    
    epsilon = 10^(-15)
    %approximate solution at different values of n
    for n = [5,9,13]
        disp('==============================================');
        fprintf('At n = %d\n',n);
        %Build varndermonde matrix A
        A = ones(m,n);
        for i = 2:n
            A(:,i) = A(:,i-1).*x;
        end
        
        %The exact solution c
        c = zeros(n,1);
        c(1:2) = 1;
        
        %The approximate solution
        [U,S,V] = svd(A);
        V = V';
        disp('After doing SVD for A, we have the system S*gamma=beta. where beta = U''y and gamma = V''c.');
        beta = U'*y
        beta_bar = beta;
        
        %this section is responsible for removing the small eigenvectors
        %from beta. removed as per request of the professor.
        %w = length(beta);
        %rel = norm(beta_bar-beta)/norm(beta);
        %while rel < epsilon
        %    temp = beta(w);
        %    beta_bar(w) = 0;
        %    rel = norm(beta_bar-beta)/norm(beta);
        %    w = w -1;
        %end
        %beta_bar(w+1) = temp
        
        gamma = S\beta_bar;
        c_hat = V'*gamma
        
        fprintf('The relative error is %e\n',norm(c_hat-c)/norm(c));
        fprintf('The residual error is %e\n',norm(A*c_hat-y));
        disp('The error in this case is much smaller than that in normal equation method.');
    end
    diary off
end