function run_class_recognition_double_split(class_dim,split_dim,pd_range,classification,generate_new_learn_model,generate_new_split)
%LoadDataAVLetterFiles({'*'},{'*'},{'*'},0);
%clear;
global database
disp '==========================================='
load('configuration');
load('all_GCF');
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
%     if i==1
%         continue;
%     end
    test_values = database.dimensions(d).values(i);
    test_GCF = query_tensor(GCF,{split_dim},{test_values});
    
    sz = size(test_GCF);
    %sz = ;
    
    lbls = get_unfolding_label(c,sz(1:end-1));
    test_GCF = reshape(test_GCF , [],sz(end));
    singleRun(c,split_dim,test_values,test_GCF ,lbls,pd_range,classifier,generate_new_learn_model,generate_new_split);
    disp '-----------------------------------'
end


for i = pd_range
    get_average_files(sprintf('Test_KPLS_*_KplsDim-%d.mat',i));
end
% for v = database.dimensions(d).values'
%     get_average_files(sprintf('Test_KPLS_*_%s-%s_*.mat',split_dim,v{1}));
% end
end


function singleRun(class_dim,split_dim_name,testValues,test_Sector ,sector_lbls,pd_range,classifier,generate_new_learn_model,generate_new_split)

run_stamp = sprintf('%s-%s',split_dim_name,cell2string(testValues));
fid = fopen('results.txt','a');
fprintf(fid,'\n==========================');
fprintf(fid,'\n single run id = %s',run_stamp);

mm = size(test_Sector,1);

    

for pd = pd_range
    %class_dim = 2;
    disp 'try dimension'
    pd
    recog_rate=0;
    recog_rate_unit=zeros(length(unique(sector_lbls)),1);
    unit_conflict=[];
    %run_stamp1 = sprintf('DIM-%d_%s_NEWDIST_KplsDim-%d',class_dim,run_stamp,pd);
    run_stamp1 = sprintf('DIM-%d_%s_KplsDim-%d',class_dim,run_stamp,pd);
%     data_file = sprintf('Test_KPLS_%s',run_stamp1);
%     if exist([data_file '.mat'],'file') && ~strcmp(classifier,'SVM')
%         continue;
%     end
        
    for ii = 1:mm
        f = ones(1,mm);
        f(ii) = 0;
        f = logical(f);

        train_GCF = test_Sector(f,:);
        train_labels = sector_lbls(f);
        test_GCF = test_Sector(ii,:);
        test_labels = sector_lbls(ii);
        splittedDataFile = 'train_test_data';
        save(splittedDataFile,'train_GCF' ,'train_labels','test_GCF' ,'test_labels' );
        
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
        [rr,rru,uc] = Test_KPLS(run_stamp1,model_file,splittedDataFile,class_dim,classifier,0);
        recog_rate = recog_rate + rr;
        
        recog_rate_unit(test_labels) = recog_rate_unit(test_labels) + rr;
%         if isempty(recog_rate_unit)
%             recog_rate_unit = rr;
%         else
%             recog_rate_unit = recog_rate_unit + rr;
%         end
%         
%         if isempty(unit_conflict)
%             unit_conflict = rr;
%         else
%             unit_conflict = unit_conflict + rr;
%         end
        
        
    end
    recog_rate = recog_rate/mm;
    %recog_rate_unit = recog_rate_unit/mm;
    save(sprintf('%s_%s','Test_KPLS', run_stamp1),'recog_rate','recog_rate_unit','unit_conflict');
    
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