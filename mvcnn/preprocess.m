load('data/split_GCF_split-test2.mat')
train_data = cell2mat({train_GCF.C})';
train_labels = cell2mat({train_GCF.class})';
test_data = cell2mat({test_GCF.C})';
test_labels = cell2mat({test_GCF.class})';
[coeff2,score2,latent2] = pca(train_data);