load('estimated_se18');%'estylevectors','eexpvectors','expression_indexes');

[totalDataNum,classNum]  = size(expression_indexes);

%%% collect feature vectors and corresponding class index %%%%%%%
featureDim = size(eexpvectors,4);
PERSON_EXP_NUM = 5;
featurevectors = zeros(totalDataNum*PERSON_EXP_NUM, featureDim);
trainFeatures = [];
trainClassIndexes = [];
trainNum = totalDataNum-1;
for i = 1:trainNum
    current_exp_indexes = expression_indexes(i,:);
    current_exp_feature_vector = squeeze(eexpvectors(i,:,:));
    current_exp_class_index = 1:6;
    null_index = find(current_exp_indexes<=0)
    current_exp_feature_vector(null_index,:) = [];
    current_exp_class_index(null_index) = [];
    trainFeatures = [trainFeatures; current_exp_feature_vector];
    trainClassIndexes = [trainClassIndexes current_exp_class_index];
end


testFeatures = [];
testClassIndexes = [];
for i = trainNum+1:totalDataNum
    current_exp_indexes = expression_indexes(i,:);
    current_exp_feature_vector = squeeze(eexpvectors(i,:,:));
    current_exp_class_index = 1:6;
    null_index = find(current_exp_indexes<=0);
    current_exp_feature_vector(null_index,:) = [];
    current_exp_class_index(null_index) = [];
    testFeatures = [testFeatures; current_exp_feature_vector];
    testClassIndexes = [testClassIndexes current_exp_class_index];
end

   
%%%%%%%%%%%%% SVM training %%%%%%%%%%%%%%%%%%%%
A = dataset(trainFeatures,trainClassIndexes');
% [m,k], a set of m datavectors of k features
% [m,n]  labels for each of the datavectors either in string or in numbers


% ---- support vector machine
p=0.25; % deviation in gaussian radial basis function
c= 1.0; % overlap of classifier
[W,J] = svc(A,'r',p,c)

% % ---- support vector machine
% p=1; % deviation in gaussian radial basis function
% c= 0.8; % overlap of classifier
% [W,J] = svc(A,'d',p)


% trainSet = dataset(trainFeatures);
% sb = W*trainSet;
% sc = classd(sb);
% 
% [index] = find (sc == trainClassIndexes');
% length(index)


testSet = dataset(testFeatures);
sb = W*testSet;
sc = classd(sb);

[index] = find (sc == testClassIndexes');
length(index)





scis = [];
%for i=1:7
    classfreq = hist(reshape(sc, dataSetGroup, classNum),classNum)
    %classfreq = hist(reshape(sc((i-1)*dataSetNum+1:i*dataSetNum,1), dataSetGroup, classNum),classNum)
    [cn ci]=max(classfreq)
    scis = [scis ; ci];
    %end

% p
% C
scis
% 







% % ---- k nearest neighbour classifier
% [W,knn,e,ek] = knnc(A)
% 
% index = 2;
% %testSet = dataset(contentVectors((index-1)*dataSetNum+1:index*dataSetNum,:),classIndex);
% testSet = dataset(contentVectors(1*dataSetNum+1:8*dataSetNum,:),repmat(classIndex,7,1));
% b = W*testSet;
% c = classd(b);
% 
% cis = [];
% for i=1:7
%     classfreq = hist(reshape(c((i-1)*dataSetNum+1:i*dataSetNum,1), dataSetGroup, classNum),classNum)
%     [cn ci]=max(classfreq)
%     cis = [cis ; ci];
% end
% 
% cis
% 

% 
% % ---- support vector machine
% p=10.0;
% C=0.8;
% % [W,J] = svc(A,'r',p,C);
% [W,J] = svc(A,'r',p);


return



% k nearest neighborhood classifier test
%load trial_xx_rec2_d40_I20.mat
%load ../gaitdata/trial_xx_rec2_d30_I20.mat
%load ../gaitdata/trial_xx_rec2_d50_I20.mat
% load ../gaitdata/trial_xx_rec2_d60_I50.mat
%load ../gaitdata/trial_xx_rec2_d100_I50.mat

load ../gaitdata50/trial38_xx_sym_dfull_GAR_I1.mat % trainVectors
load ../gaitdata50/trial38_xx_sym_dfull_CBL_I1.mat  % contentVectors


[totalDataNum dataDim] = size(contentVectors);
dataSetGroup = 7;
classNum = 38;
dataSetNum = totalDataNum/dataSetGroup;
sampleCycle = dataSetNum/classNum;

% classIndex = reshape(repmat(result(1,:),sampleCycle,1),1,dataSetNum)';
tempIndex = 1:38;
classIndex = reshape(repmat(tempIndex,dataSetGroup,1),1,dataSetGroup*classNum)';

% A = dataset(trainVectors,classIndex);
% 
% % % ---- k nearest neighbour classifier
% % [W,knn,e,ek] = knnc(A)
% % 
% % index = 2;
% % %testSet = dataset(contentVectors((index-1)*dataSetNum+1:index*dataSetNum,:),classIndex);
% % testSet = dataset(contentVectors(1*dataSetNum+1:8*dataSetNum,:),repmat(classIndex,7,1));
% % b = W*testSet;
% % c = classd(b);
% % 
% % cis = [];
% % for i=1:7
% %     classfreq = hist(reshape(c((i-1)*dataSetNum+1:i*dataSetNum,1), dataSetGroup, classNum),classNum)
% %     [cn ci]=max(classfreq)
% %     cis = [cis ; ci];
% % end
% % 
% % cis
% % 
% 
% 
% % ---- support vector machine
% p=10.0;
% C=0.8;
% % [W,J] = svc(A,'r',p,C);
% [W,J] = svc(A,'r',p);

load ../gaitdata50/SVM_weight38_j266; %W','J','trainVectors'



testSet = dataset(contentVectors,classIndex);
%testSet = dataset(contentVectors(1*dataSetNum+1:8*dataSetNum,:),repmat(classIndex,7,1));
sb = W*testSet;
sc = classd(sb);

scis = [];
%for i=1:7
    classfreq = hist(reshape(sc, dataSetGroup, classNum),classNum)
    %classfreq = hist(reshape(sc((i-1)*dataSetNum+1:i*dataSetNum,1), dataSetGroup, classNum),classNum)
    [cn ci]=max(classfreq)
    scis = [scis ; ci];
    %end

% p
% C
scis
% 



load('TensorEXP_p18e6-5_cf');%,'GCFM','rbf_cent','CENTER_NUM','DATA_DIM','subject_indexes','expression_indexes');


save('estimated_se18','estylevectors','eexpvectors','expression_indexes');