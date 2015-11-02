function GCF = getCoeffTensor_normalized(Seqs,ncent, save_to_file)
    %CF=cell(size(Seqs));
    % for i=1:size(trainSeqs,1)
    %     for j=1:size(trainSeqs,2)
    %             % embed on a circle
    %             frames=trainSizes(i,j);%angletrain{seq};
    %             t=[1:frames]'*2*pi/frames;
    %             P=[cos(t) sin(t)];
    %             % learn mapping between the manifolds and the input space
    %             CF{i,j}=learnmapping_grbf(trainSeqs{i,j}',P);
    %             %Y = Si X, where Y is the input sequence (trainSeqs), X is the embedded
    %             %space (P) and Si is the mapping between both spaces.
    %     end
    % end
    debug = 1;
    if exist('ncent','var') && ncent>0
        icent=((1:ncent)-1)'*2*pi/ncent;
        cent = [cos(icent) sin(icent)]; 
    else
        cent = [];
    end
    
    sz = size(Seqs);
    l = learnmapping_grbf_regulize(Seqs{1,1,1}',cent);
    GCF = zeros([sz,size(l)]);
    
    for i=1:size(Seqs,1)
        for j=1:size(Seqs,2)
            for k = 1:size(Seqs,3)
    %             % embed on a circle
    %             frames=trainSizes(i,j,k);%angletrain{seq};
    %             t=[1:frames]'*2*pi/frames;
    %             P=[cos(t) sin(t)];
                % learn mapping between the manifolds and the input space
                %CF{i,j,k}=learnmapping_grbf_regulize(Seqs{i,j,k}');
                GCF(i,j,k,:,:)=learnmapping_grbf_regulize(Seqs{i,j,k}',cent);
                %Y = Si X, where Y is the input sequence (trainSeqs), X is the embedded
                %space (P) and Si is the mapping between both spaces.
            end
        end
    end
    
    %clear trainSeqs;
    clear l;
    clear Seqs;
    %% normalization
    %mean
    GCF = reshape(GCF,prod(sz),[]);
    l = mean(GCF);
    GCF = bsxfun(@minus, GCF, l);
    if debug
        h = figure; colormap(gray)
        imagesc(cosine_similarity(GCF,GCF));
        print(h, '-dpng', sprintf('ncents_%d_GCF_mean',ncent));
    end
    %variance
    v = sqrt(var(GCF));
    GCF =  bsxfun(@rdivide,GCF,v);
    %remove NaN elements from GCF, which is the columns with variance zero
    [~,y] = find(isnan(GCF));
    y = unique(y);
    ind = true(1,size(GCF,2));
    ind(y) = false;
    GCF = GCF(:,ind);
    
    
    %% normalize rows' norm
    nrm = sqrt(diag(GCF*GCF'));
    nrm = nrm(:);
    GCF = bsxfun(@rdivide, GCF,nrm);
    
    %imagesc(cosine_similarity(GCF,GCF));
    %print(h, '-deps', sprintf('ncents_%d_GCF_mean_norm',ncents));
    %% final reshaping
    GCF = reshape(GCF,sz(1),sz(2),sz(3),[]);
    
    if exist('save_to_file','var') && save_to_file
        filename = sprintf('all_GCF_ncents_%d',ncent);
        save(filename,'GCF');
    end
    % %GCF = cell2mat(CF);
    % GCF = zeros([size(CF),size(CF{1,1})]);
    % for i=1:size(GCF,1)
    %     for j=1:size(GCF,2)
    %         
    %             GCF(i,j,:,:) = CF{i,j};
    %         
    %     end
    % end
    
%     for i=1:size(GCF,1)
%         for j=1:size(GCF,2)
%             for k = 1:size(GCF,3)
%                 GCF(i,j,k,:,:) = CF{i,j,k};
%             end
%         end
%     end
%    clear CF;
end