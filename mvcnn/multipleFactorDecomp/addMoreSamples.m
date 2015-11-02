function B = addMoreSamples(A,dim,n)
    sigma = 0.005;
    mu = 0;
    A2d = unfold(A,dim);
    os = size(A2d,1);
    C = zeros(n*os,size(A2d,2));
    C(1:os,:) = A2d;
    for i = 1:n-1
        S = normrnd(mu,sigma,size(A2d));
        C((i*os+[1:os]),:) = S;
    end
    ns = size(A);
    ns(dim) = ns(dim)*n;
    B = fold(C,dim,ns);
end