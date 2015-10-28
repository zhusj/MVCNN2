load('./data/split_GCF_split-test.mat')
train_data = cell2mat({train_GCF.C})';
% train_labels = cell2mat({train_GCF.class})';
% test_data = cell2mat({test_GCF.C})';
% test_labels = cell2mat({test_GCF.class})';
% [coeff2,score2,latent2] = pca(train_data);
% [coeff2,score2,latent2] = princomp(train_data);
% [residuals,reconstructed] = pcares(test_data,4096);

% 
d = 500;
train_mean = mean(train_data);
Basis_ = bsxfun(@minus,train_data,train_mean);
[Basis_, S, V] = myPCA(Basis_,d);
s = diag(S);
W_d = V'*diag(1./s);% this is the projection matrix which I believe can be used for initializing FC7