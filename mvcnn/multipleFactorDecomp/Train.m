function storeDataFile = Train(run_stamp,p_load_train_data,repeats,letters,subjects)
%clear
%%preprocessing
diary 'train_log.txt'
echo on
%%load data
%[trainSeqs,trainSizes,trainLetters,trainSubjects] = LoadAVLetterTrainingData;

if ~exist('repeats','var')
    repeats = {'1','2'}
end

if ~exist('letters','var')
    letters = {'*'}
end

if ~exist('subjects','var')
    subjects = {'*'}
end

%trainDataFile = sprintf('load_TRAIN_data_%s',run_stamp);
GCF = LoadDataAVLetterFiles(letters,subjects,repeats,p_load_train_data);


% if exist('p_load_train_data','var') && p_load_train_data
%     load(trainDataFile);
% else
%     [tSeqs,~,tLetters,tSubjects] = LoadDataAVLetterFiles(letters,subjects,repeats,trainDataFile);
% end




% in the sampledata.mat, there are 5 sequences (categories) for training
% and testing

%%processing
x={};%embedding coordinates
w={};%weights of each class style
a={};%estimated viewpoint
er={};%angle error


% learn mapping between the manifolds and the input space
%disp('learn mapping between the manifolds and the input space');
%GCF = getCoeffTensor(tSeqs);
%clear trainSeqs


storeDataFile = sprintf('AVLetter_TRAIN_GCF_%s',run_stamp);
save(storeDataFile,'GCF');%,'tLetters','tSubjects');

echo off
diary off
return;
%% at this point construction of the data and coefficient data is finished

% trainLetters = reshape(repmat(trainLetters,1,size(trainSeqs,2)),[],1);
% trainSubjects = reshape(repmat(trainSubjects',size(trainSeqs,2),1),[],1);

% trainSeqs = reshape(trainSeqs,size(trainSeqs,1)*size(trainSeqs,2),[]);
% trainSizes = reshape(trainSizes,size(trainSizes,1)*size(trainSizes,2),[]);
% trainSeqs = reshape(trainSeqs,size(trainSeqs,1),[]);
% trainSizes = reshape(trainSizes,size(trainSizes,1),[]);

%load('AVLetter_D1-Letters_D2-Variations_D3-Subjects_D4-5-SequenceFeatures');
load('AVLetter_D1-Letters_D2-reptns(1-2)_D3-Subjects_D4-5-SequenceFeatures');
% 
% 
% trainLetters = reshape(repmat(trainLetters,1,size(GCF,2)),[],1);
% GCF = reshape(GCF,size(GCF,1)*size(GCF,2),size(GCF,3),size(GCF,4),[]);
% B = unfold(GCF,1);

trainSubjects = reshape(repmat(trainSubjects',size(GCF,2),1),[],1);
GCF = reshape(GCF,size(GCF,1),size(GCF,2)*size(GCF,3),size(GCF,4),[]);
B = unfold(GCF,2);
clear GCF;

%decomposititon
%style vectors are the rows of V
[U,S,V]=svd(B',0);
%[V,S] = eig(B);
%V1 =  compute_mapping(V, 'Isomap',3);
V1 = V(:,1:3);
%V1 = U(:,1:3);
%plotting(V1,trainLetters);
plotting(V1,trainSubjects);
 
% solving for style and content of test dataset
disp('solving for style and content');

%%Testingsize(U)
%[testSeqs,testSizes,testLabel,testSubject,TNseq] = LoadAVLetterTestData
% for nseq=1:1:TNseq
%     
%     
%     message=sprintf('solvig for style and content of instance %d',nseq);
%     disp(message); 
%     
%     nff=size(testSeqs{nseq},2);
%     
%     for nf=1:1:nff
%         y=testSeqs{nseq}(:,nf);
%         [x{nseq}(:,nf),w{nseq}(:,nf),a{nseq}(:,nf)]=solv4sc(U,S,V,y,8*8*9,Nb+d+1,cent);
% %         eee=a{nseq}(:,nf)*180/pi-angletest{nseq}(nf);
% %         if eee>180
% %             eee=eee-360;
% %         end
% %         if eee<-180
% %             eee=eee+360;              
% %         end
% %        er{nseq}(nf)=eee;
%     end
%         
% end
disp('The End')