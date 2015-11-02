function GCF = getCoeffTensor_normalized(Seqs, save_to_file,varargin)

    debug = 1;
    [GCF,ncents] = getCoeffTensor(Seqs,save_to_file,varargin{:});
    
    normalized = true;
    GCF_sz = size(GCF);
    %clear trainSeqs;
    clear Seqs;
    %% normalization
    %mean
    %GCF_sz = [GCF_sz(1:end-2),prod(l_sz)];
    GCF = reshape(GCF,[],GCF_sz(end));
    l = mean(GCF);
    GCF = bsxfun(@minus, GCF, l);
    if debug
        h = figure; colormap(gray)
        imagesc(cosine_similarity(GCF,GCF));
        print(h, '-dpng', sprintf('ncents_%d_GCF_mean',ncents));
    end
    %variance
    size_bfr_varinace = size(GCF);
    v = sqrt(var(GCF));
    GCF =  bsxfun(@rdivide,GCF,v);
    %remove NaN elements from GCF, which is the columns with variance zero
    [~,y] = find(isnan(GCF));
    y = unique(y);
    ind = true(1,size(GCF,2));
    ind(y) = false;
    GCF = GCF(:,ind);
    size_afr_varinace = size(GCF);
    
    %% normalize rows' norm
    nrm = sqrt(diag(GCF*GCF'));
    nrm = nrm(:);
    GCF = bsxfun(@rdivide, GCF,nrm);
%     
    %imagesc(cosine_similarity(GCF,GCF));
    %print(h, '-deps', sprintf('ncents_%d_GCF_mean_norm',ncents));
    GCF_sz(end) = size(GCF,2);
    %% final reshaping
    GCF = reshape(GCF,GCF_sz);
    
    if exist('save_to_file','var') && save_to_file
%         %filename = sprintf('all_GCF_ncents_%d',ncents);
        filename = sprintf('all_GCF');
        save(filename,'GCF','normalized','size_bfr_varinace','size_afr_varinace','-append');
%         s = whos('GCF');
%         
%         if s.bytes > 2*1024^3
%             disp 'the matrix GCF is huge'
%             save(filename,'GCF','ncents','normalized','size_bfr_varinace','size_afr_varinace','-v7.3');
%         else
%             save(filename,'GCF','ncents','normalized','size_bfr_varinace','size_afr_varinace');
%         end
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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