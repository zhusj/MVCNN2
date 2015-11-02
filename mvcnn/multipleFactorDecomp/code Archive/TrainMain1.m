function GCF = TrainMain2
clear
%%preprocessing
diary 'train_log.txt'
echo on
%%load data
%[trainSeqs,trainSizes,trainLetters,trainSubjects] = LoadAVLetterTrainingData;
subjects = {'*'};
letters = {'*'};
repeats = {'2','3'};
[trainSeqs,trainSizes,trainLetters,trainSubjects] = LoadDataAVLetterFiles(letters,subjects,repeats);



% in the sampledata.mat, there are 5 sequences (categories) for training
% and testing

%%processing
x={};%embedding coordinates
w={};%weights of each class style
a={};%estimated viewpoint
er={};%angle error


% % embed on a circle
% for seq=1:length(trainSeqs)
%     frames=trainSizes(seq);%angletrain{seq};
%     ti = 1:frames;
%     t=ti*2*pi/frames;
%     P{seq}=[cos(t') sin(t')];
% end

% learn mapping between the manifolds and the input space
disp('learn mapping between the manifolds and the input space');
% M=10; % mapping centers, Need to know if this number affects the performance???????
% ti1=[1:M];
% t1=ti1*2*pi/M;
% cent=[cos(t1') sin(t1')];
% 
% clear('ti','t1','ti1');
% 
% Nb=length(cent);
% d=2;




%CF=cell(length(trainLetters),length(trainSubjects));
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
for i=1:size(trainSeqs,1)
    for j=1:size(trainSeqs,2)
        for k = 1:size(trainSeqs,3)
%             % embed on a circle
%             frames=trainSizes(i,j,k);%angletrain{seq};
%             t=[1:frames]'*2*pi/frames;
%             P=[cos(t) sin(t)];
            % learn mapping between the manifolds and the input space
            CF{i,j,k}=learnmapping_grbf_regulize(trainSeqs{i,j,k}');
            %Y = Si X, where Y is the input sequence (trainSeqs), X is the embedded
            %space (P) and Si is the mapping between both spaces.
        end
    end
end
clear trainSeqs

% %GCF = cell2mat(CF);
% GCF = zeros([size(CF),size(CF{1,1})]);
% for i=1:size(GCF,1)
%     for j=1:size(GCF,2)
%         
%             GCF(i,j,:,:) = CF{i,j};
%         
%     end
% end
GCF = zeros([size(CF),size(CF{1,1,1})]);
for i=1:size(GCF,1)
    for j=1:size(GCF,2)
        for k = 1:size(GCF,3)
            GCF(i,j,k,:,:) = CF{i,j,k};
        end
    end
end
clear CF;
%save('AVLetter_D1-Letters-Variations_D2-Subjects_D3-4-SequenceFeatures','GCF','trainLetters','trainSubjects');
%save('AVLetter_D1-Letters_D2-Variations-Subjects_D3-4-SequenceFeatures','GCF','trainLetters','trainSubjects');
%save('AVLetter_D1-Letters_D2-reptns(1-2)_D3-Subjects_D4-5-SequenceFeatures','GCF','trainLetters','trainSubjects');
save('AVLetter_D1-Letters_D2-reptns(2-3)_D3-Subjects_D4-5-SequenceFeatures','GCF','trainLetters','trainSubjects');

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