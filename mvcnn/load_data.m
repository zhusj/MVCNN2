clear
data_path = '/media/DATA/mvcnn/features/';
model = 'ModelNet40v1-imagenet-vgg-m-finetuned-ModelNet40v1-BS60_AUGnone-finetuned-ModelNet40v1-fine_tuned_vgg_add_12_branch_2048_after_fc6_LRw4_ave_pooling_leaky0.01_8epochs_rmfc7_MVconv5-none';

path = strcat(data_path, model, '/NORM0/');

files = dir(fullfile(path,'fc_b*.mat'));

for i = 1:length(files)
    fileLists{i} = files(i).name;
end

[S,INDEX] = sort_nat(fileLists);

for i = 1:length(S)
    data{i} = load(fullfile(path,S{i}),'x');
end

for i=1:12
    x(:,:,i) = data{i}.x;
end

totalDist = 0;
for i =1:3983
    x1(:,:) = x(i,:,:);
    x1 = x1/norm(x1(:));
    dist = pdist(x1');
    totalDist = totalDist + dist;
end

sum(totalDist)