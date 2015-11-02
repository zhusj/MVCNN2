function run
clear;
disp '-----------------------------------'
trainRep = {'1','2'}
testRep = {'3'}
singleRun(trainRep,testRep,1);
%%%%%%%%%%%%%%%
disp '-----------------------------------'
clear;
trainRep = {'1','3'}
testRep = {'2'}
singleRun(trainRep,testRep,0)
%%%%%%%%%%%%%%%%%%%%%%%
disp '-----------------------------------'
clear;
trainRep = {'2','3'}
testRep = {'1'}
singleRun(trainRep,testRep,0)
end


function singleRun(trainRep,testRep,load_saved_data)
if ~exist('load_saved_data','var')
    load_saved_data = 0;
end
disp '**TRAIN'
trainData = Train(trainRep,load_saved_data);
%trainData = sprintf('AVLetter_D1-Letters_D2-reptns_%s_D3-Subjects_D4-5-SequenceFeatures',cell2string(trainRep));
%load(sprintf('AVLetter_D1-Letters_D2-reptns_%s_D3-Subjects_D4-5-SequenceFeatures',cell2string(repeats)));
disp '**LEARN'
model_file = learn_HOSVD(trainData,trainRep,1);
%model_file = sprintf('AVLetter_Learn_HOSVD_52-10-(%s)',cell2string(trainRep));
%model_file = learn_LDA(trainData,trainRep,1);
%model_file = sprintf('AVLetter_Learn_LDA_52-10-(%s)',cell2string(trainRep));

disp '**TEST'
Test_knn_HOSVD(model_file,testRep,load_saved_data)
%Test_knn_LDA(model_file,testRep,load_saved_data)
%plot_letter_confusion_matrix('results_test_NEW');
end