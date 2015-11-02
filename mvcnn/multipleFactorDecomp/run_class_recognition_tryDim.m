function run_class_recognition_tryDim(class_dim,split_dim,pd_range,classification,generate_new_learn_model,generate_new_split)
%LoadDataAVLetterFiles({'*'},{'*'},{'*'},0);
%clear;
global database
disp '==========================================='
load('configuration');

if nargin < 1
    split_dim = 'subjects';
    class_dim = 'words';
end
if ~exist('generate_new_learn_model','var')
    generate_new_learn_model = 0;
end

if ~exist('generate_new_split','var')
    generate_new_split = 0;
end

if ~exist('classification','var')
    classifier='';
else
    classifier = classification;
end

if ~exist('pd_range','var')
    pd_range = 25%40:10:130;%[15:5:50 60:10:100];
end

% if strcmp(classifier,'SVM')
%     classifier = ['_' classifier];
% end

d = strmatch(split_dim, {database.dimensions.name}, 'exact');
c = strmatch(class_dim, {database.dimensions.name}, 'exact');
for i = 1:length(database.dimensions(d).values)
    f = ones(1,length(database.dimensions(d).values));
    f(i) = 0;
    trainValues = database.dimensions(d).values(logical(f));
    testValues = database.dimensions(d).values(i);
    singleRun(c,split_dim,trainValues,testValues,pd_range,classifier,generate_new_learn_model,generate_new_split);
    disp '-----------------------------------'
end


for i = pd_range
    get_average_files(sprintf('Test_KPLS_*_KplsDim-%d.mat',i));
end
end


function singleRun(class_dim,split_dim_name,trainValues,testValues,pd_range,classifier,generate_new_learn_model,generate_new_split)

run_stamp = sprintf('%s-%s',split_dim_name,cell2string(testValues));
fid = fopen('results.txt','a');
fprintf(fid,'\n==========================');
fprintf(fid,'\n single run id = %s',run_stamp);

disp '**TRAIN'
splittedDataFile = sprintf('split_GCF_%s',run_stamp);
if ~exist([splittedDataFile '.mat'],'file') || generate_new_split
    splittedDataFile = split_data(run_stamp,split_dim_name,trainValues,testValues);
end

sprintf('**LEARN:');
%model_file = learn_HOSVD(run_stamp,trainData,0);
%model_file = learn_HOSVD_KPLS(run_stamp,trainData,1);
%model_file = sprintf('AVLetter_Learn_HOSVD_%s',run_stamp);
%model_file = learn_LDA(trainData,trainRep,1);
%model_file = sprintf('AVLetter_Learn_LDA_52-10-(%s)',cell2string(trainRep));

for pd = pd_range
    %class_dim = 2;
    disp 'try dimension'
    pd
    %run_stamp1 = sprintf('DIM-%d_%s_NEWDIST_KplsDim-%d',class_dim,run_stamp,pd);
    run_stamp1 = sprintf('DIM-%d_%s_KplsDim-%d',class_dim,run_stamp,pd);
    sprintf('**LEARN:');
    %model_file = learn_KPLS(run_stamp1,trainData,class_dim,0,pd);
    
    
    %model_file = sprintf('Learn_KPLS_SVM_%s',run_stamp1);
    model_file = sprintf('learn_KPLS_SVM_%s',run_stamp1);
    
    if ~exist([model_file '.mat'],'file') || generate_new_learn_model
        model_file = learn_KPLS_SVM(run_stamp1,splittedDataFile,class_dim,pd,0,split_dim_name);
    end
    

    
    %model_file = sprintf('AVLetter_Learn_KPLS_%s',run_stamp);

    disp '**TEST'
    %Test_knn_HOSVD(run_stamp,model_file,load_saved_data,testRep,{'*'},{'*'})
    %Test_knn_LDA(model_file,testRep,load_saved_data)
    %plot_letter_confusion_matrix('results_test_NEW');
%     if strcmp(classifier,'_SVM')
%         recog_rate = Test_KPLS_SVM(run_stamp1,model_file,splittedDataFile,class_dim);
%     else
%         recog_rate = Test_KPLS(run_stamp1,model_file,splittedDataFile,class_dim);
%     end
    recog_rate = Test_KPLS(run_stamp1,model_file,splittedDataFile,class_dim,classifier);
    fprintf(fid,'\n KPLS_dimension(%d)=> recognition rate %f',pd,recog_rate);
    load(model_file,'pls_dim');
    if pls_dim < pd
        break;
    end
    %delete(model_file);
end
fprintf(fid,'\n==== end of single run===');

fclose('all');
end