function run_class_recognition(class_dim,split_dim,generate_new_learn_model)
%LoadDataAVLetterFiles({'*'},{'*'},{'*'},0);
%clear;
disp '==========================================='
load('configuration');
if nargin < 1
    split_dim = 'subjects';
    class_dim = 'words';
end

if ~exist('generate_new_learn_model','var')
    generate_new_learn_model = 0;
end

d = strmatch(split_dim, {database.dimensions.name}, 'exact');
c = strmatch(class_dim, {database.dimensions.name}, 'exact');
for i = 1:length(database.dimensions(d).values)
    f = ones(1,length(database.dimensions(d).values));
    f(i) = 0;
    trainValues = database.dimensions(d).values(logical(f));
    testValues = database.dimensions(d).values(i);
    singleRun(c,split_dim,trainValues,testValues,generate_new_learn_model);
    disp '-----------------------------------'
end

get_average_files(sprintf('test_results_DIM-%d*.mat',c));

end


function singleRun(class_dim,split_dim_name,trainValues,testValues,generate_new_learn_model)

run_stamp = sprintf('%s_%s',split_dim_name,cell2string(testValues));
fid = fopen('results.txt','a');
fprintf(fid,'\n==========================');
fprintf(fid,'\n single run id = %s',run_stamp);

disp '**TRAIN'
splittedDataFile = sprintf('split_GCF_%s',run_stamp);
if ~exist([splittedDataFile '.mat'],'file')
    splittedDataFile = split_data(run_stamp,split_dim_name,trainValues,testValues);
end


sprintf('**LEARN:');
%model_file = learn_HOSVD(run_stamp,trainData,0);
%model_file = learn_HOSVD_KPLS(run_stamp,trainData,1);
%model_file = sprintf('AVLetter_Learn_HOSVD_%s',run_stamp);
%model_file = learn_LDA(trainData,trainRep,1);
%model_file = sprintf('AVLetter_Learn_LDA_52-10-(%s)',cell2string(trainRep));


%class_dim = 2;
run_stamp1 = sprintf('DIM-%d_%s',class_dim,run_stamp);
sprintf('**LEARN:');

model_file = sprintf('Learn_KPLS_SVM_%s',run_stamp1);
if ~exist([model_file '.mat'],'file') || generate_new_learn_model
    model_file = learn_KPLS_SVM(run_stamp1,splittedDataFile,class_dim,0);
end

disp '**TEST'
%Test_knn_HOSVD(run_stamp,model_file,load_saved_data,testRep,{'*'},{'*'})
%Test_knn_LDA(model_file,testRep,load_saved_data)
%plot_letter_confusion_matrix('results_test_NEW');

recog_rate = Test_KPLS_SVM(run_stamp1,model_file,splittedDataFile,class_dim);
fprintf(fid,'\n For %s: recognition rate %f',run_stamp1,recog_rate);
delete(model_file);

fprintf(fid,'\n==== end of single run===');

fclose('all');
end