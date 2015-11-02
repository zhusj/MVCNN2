function GCF = getCoeffTensor(Seqs)
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
    
    l = learnmapping_grbf_regulize(Seqs{1,1,1}');
    GCF = zeros([size(Seqs),size(l)]);
    
    for i=1:size(Seqs,1)
        for j=1:size(Seqs,2)
            for k = 1:size(Seqs,3)
    %             % embed on a circle
    %             frames=trainSizes(i,j,k);%angletrain{seq};
    %             t=[1:frames]'*2*pi/frames;
    %             P=[cos(t) sin(t)];
                % learn mapping between the manifolds and the input space
                %CF{i,j,k}=learnmapping_grbf_regulize(Seqs{i,j,k}');
                GCF(i,j,k,:,:)=learnmapping_grbf_regulize(Seqs{i,j,k}');
                %Y = Si X, where Y is the input sequence (trainSeqs), X is the embedded
                %space (P) and Si is the mapping between both spaces.
            end
        end
    end
    %clear trainSeqs;
    clear l;
    GCF = reshape(GCF,size(GCF,1),size(GCF,2),size(GCF,3),[]);
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