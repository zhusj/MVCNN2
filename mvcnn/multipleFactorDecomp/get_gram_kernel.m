function K =  get_gram_kernel(A,B)
    %K = A*B';
    %K = exp(-0.5*K);
    %K = exp(-0.05*dist2(A,B));
    %K = exp(-0.05*pdist2(A,B));
    %K = exp(-1*pdist2(A,B,'cosine'));
    % K = bsxfun(@minus,K, mean(K));
    % K = bsxfun(@rdivide,K, sqrt(var((K))));
    
%     if ~exist('B','var')
%         
%         K = dist2(A,A);
%     else
%         K = dist2(A,B);
%     end
%     mm = mode(abs(K(:)));
%     K = exp(-K/mm);
%     return;
    %% cosine similarity
    if ~exist('B','var')
        K = cosine_similarity(A);
    else
        K = cosine_similarity(A,B);
    end
    return;
    %% Grassmannian Kernel
    load('all_GCF','ncents');
    sz = size(A);
    C1 = zeros(sz(end)/ncents,ncents,sz(1));
    for i = 1:sz(1)
        C1(:,:,i) = grammScmidt(reshape(A(i,:),ncents,[])');
    end
    if ~exist('B','var')
        K = Compute_Grassmann_Kernel(C1);
    else
        sz = size(B);
        C2 = zeros(sz(end)/ncents,ncents,sz(1));
        for i = 1:sz(1)
            C2(:,:,i) = grammScmidt(reshape(B(i,:),ncents,[])');
        end
        K = Compute_Grassmann_Kernel(C1,C2);
    end
    mxProj = max(max(K(:,:,1)));
    mxCC = max(max(K(:,:,2)));
    mnProj = min(min(K(:,:,1)));
    mnCC = min(min(K(:,:,2)));
    
    weight = ([mxProj, mxCC].^-1)';
    K = squeeze(tmul(K,weight,3));