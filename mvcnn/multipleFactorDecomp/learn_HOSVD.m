function storeFile = learn_HOSVD(run_stamp,filename,plot_flg)
% FILE NAME: learn_hosvd_gait_bat.m
% tensor analysis for gcf not gscf

diary 'run_log_learn_hosvd.txt'
echo on
%clear
%pack

if ~exist('plot_flg','var')
    plot_flg = 1;
end
if nargin<1
    filename = 'AVLetter_D1-Letters_D2-reptns(1-2)_D3-Subjects_D4-5-SequenceFeatures';
    %'AVLetter_Expr10-1_Gabor-set2_D1-Letters_D2-reptns(1-2)_D3-Subjects_D4-5-SequenceFeatures';
    %'AVLetter_Gabor-set2_D1-Letters_D2-reptns(1-2)_D3-Subjects_D4-5-SequenceFeatures';
    %'AVLetter_HoG_D1-Letters_D2-reptns(1-2)_D3-Subjects_D4-5-SequenceFeatures';
    %strcat('New-AVLetter_D1-Letters_D2-reptns(2-3)_D3-Subjects_D4-5-SequenceFeatures'); %GCF    
end
load(filename); 


%% remove the repeatations dimensions
% by taking the mean over all repeatations

% first unfold the last two dimensions
% s = size(GCF);
train_seq_dim = size(GCF);
% l = length(s);
% s1 = zeros(1,l-1);
% s1(1:l-2)= s(1:l-2);
% s1(l-1) = s(l-1)*s(l);
% nGCF = reshape(GCF,s1);
% clear GCF;
% zeros([size(GCF,1),size(GCF,2),size(GCF,3),size(GCF,4)*size(GCF,5)]
% for i = 1: size(GCF,1)
%     for j = 1: size(GCF,2)
%         for k = 1: size(GCF,3)
%             nGCF(i,j,k,:) = reshape(GCF(i,j,k,:,:),1,[]);
%         end
%     end
% end

%nGCF = squeeze(mean(nGCF,2));%one option is to take the mean over all repeatations
GCF = reshape(GCF,size(GCF,1)*size(GCF,2),size(GCF,3),[]);%another option is to merge repeatation with letter to increase accuracy of letters classification
%nGCF = reshape(nGCF,size(nGCF,1),size(nGCF,2)*size(nGCF,3),size(nGCF,4));%another option is to merge repeatation with subjects
%nGCF = addMoreSamples(nGCF,2,5);

            

% rdim = [10, 3, 11, 120] % 10           3       49152 ( 12 x 4096)
rdim = size(GCF)%[10, 3, 11, 220] % 10           3       45056 ( 11 x 4096)
for d = 52
rdim(1) = d;
%rdim(end) = rdim(end)/100;
%rdim(1) = rdim(1)/2;
tic
%[c u s ] = HOSVDDR2(GCF,rdim,5);
[c u s ] = HOSVDDR2(GCF,rdim,10);
% [c u s ] = hooi_svdn(GCF,rdim);
toc
end
if plot_flg
    letter_bases = u{1};
    style_bases = u{2};
    ltr_labels = (1:train_seq_dim(1))';%[1:26,1:26]';
    ltr_labels = repmat(ltr_labels,train_seq_dim(2),1);
    style_labels = (1:train_seq_dim(3))';
    plotting(letter_bases(:,1:3),ltr_labels);
    %plotting(style_bases(:,1:3),style_labels);
end

%save('AVLetter_Train_SVD_26-10-2_reps-mean.mat','c','u','s','rdim');
%save('AVLetter_Train_SVD-merge_letters_reps.mat','c','u','s','rdim');
%save('AVLetter_HoG_Learn_SVD_26-10-(1-2)_merge-subjs-reps.mat','c','u','s','rdim');
%save('AVLetter__Gabor-set1_Learn_SVD_52-50-(1-2).mat','c','u','s','rdim');
storeFile = sprintf('AVLetter_Learn_HOSVD_%s',run_stamp);
save(storeFile,'c','u','s','rdim','train_seq_dim');
%save('AVLetter_Train_SVD_merge-reps-subjs.mat','c','u','s','rdim');
echo on
diary off
return;
%% load plot
%load('AVLetter_Train_SVD-merge_letters_reps.mat');
figure
bar(u{1})
title('Letter Vector')
figure
bar(u{2})
title('Person Style Vector')

V1 = u{1}(:,1:3);
plotting(V1,trainLetters);
V1 = u{2}(:,1:3);
plotting(V1,trainSubjects);

return;
%%% visualize tensor analysis result

V1 = u{1};

figure
plot3(V1(1,1), V1(1,2),V1(1,3),'rx');
hold on
plot3(V1(2,1), V1(2,2),V1(2,3),'go');
plot3(V1(3,1), V1(3,2),V1(3,3),'b*');
plot3(V1(4,1), V1(4,2),V1(4,3),'md');
plot3(V1(5,1), V1(5,2),V1(5,3),'k+');
plot3(V1(6,1), V1(6,2),V1(6,3),'r>');
plot3(V1(7,1), V1(7,2),V1(7,3),'g^');
plot3(V1(8,1), V1(8,2),V1(8,3),'bs');
plot3(V1(9,1), V1(9,2),V1(9,3),'m<');
plot3(V1(10,1), V1(10,2),V1(10,3),'kp');

title('Expression style component 1,2,3')
legend('person style 1','person style 2','person style 3','person style 4','person style 5','person style 6','person style 7','person style 8','person style 9','person style 10')


V2 = u{2};

% for lowindex = 1:2:4
% 	highindex = lowindex+1;
% 	figure
% 	plot(V2(1,lowindex), V2(1,highindex),'rx');
% 	hold on
% 	plot(V2(2,lowindex), V2(2,highindex),'go');
% 	plot(V2(3,lowindex), V2(3,highindex),'b*');
% 	plot(V2(4,lowindex), V2(4,highindex),'md');
% % 	plot(V2(5,lowindex), V2(5,highindex),'k+');
% 	titlest = strcat('View vector component:',int2str(lowindex),'-',int2str(highindex),' ');
%     title(titlest)
% end

figure
plot3(V2(1,1), V2(1,2),V2(1,3),'rx');
hold on
plot3(V2(2,1), V2(2,2),V2(1,3),'go');
plot3(V2(3,1), V2(3,2),V2(1,3),'b*');
% plot3(V2(4,1), V2(4,2),V2(1,3),'md');
% plot3(V2(5,1), V2(5,2),V2(1,3),'k+');
title('Expression vector component 1,2,3')


return
