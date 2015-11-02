load('all_data_files');
echo on
currentFolder = pwd;
addpath(currentFolder);
pd = [80, 90,100,130,180,200,220,250];
for n =[4,5,6,7,8,9]%4:2:10%4:2:10
    root = fullfile(pwd,sprintf('ncents_%d',n))
    mkdir(root);
    
    for rbf = {'gauss'}%,'poly'}%gaussian is always better than thinplate
        cd(root);
        folder = fullfile(pwd,sprintf('%s',rbf{1}) )
        mkdir(folder);
        cd(folder);
        addpath(pwd)
        getCoeffTensor_normalized(tSeqs,1,n,rbf{1});

        subfolder = folder%fullfile(folder,'Normal')
        %mkdir(subfolder);
        cd(subfolder)
        run_class_recognition_tryDim('words','repeats',pd,'',1,1);

        subfolder2 = fullfile(subfolder,'SVM')
        mkdir(subfolder2);
        cd(subfolder2)
        run_class_recognition_tryDim('words','repeats',pd,'SVM',0,0);

        cd(subfolder);
        delete_files_expr('split_GCF_*.mat');
        rmpath(folder);
    end
    
    cd(currentFolder);
end
echo off
    