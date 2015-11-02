function [split_data_file] = split_GCF(run_stamp,dimension_name,p_train_values,p_test_values)
%clear
%%preprocessing
diary 'train_log.txt'
echo on
load('all_GCF');
% load('configuration');
% dimension_names = {database.dimensions.name};
% indx = strmatch(dimension_name, dimension_names, 'exact');
% 
% dim_train_values = cell(length(database.dimensions),1);
% dim_test_values = cell(length(database.dimensions),1);
% for i = 1:length(database.dimensions)
%     if i == indx
%         dim_train_values{i} = p_train_values;
%         dim_test_values{i} = p_test_values;
%         
%     else
%         dim_train_values{i} = '*';
%         dim_test_values{i} = '*';
%     end
% end

if iscell(GCF) %&& (numel(size(GCF))>2)
    train_GCF = query_tensor(GCF,{dimension_name},{p_train_values});
    test_GCF = query_tensor(GCF,{dimension_name},{p_test_values});
end
if isstruct(GCF)
    l = cell2mat({GCF.(dimension_name)});
    ind_test = (ismember(l,str2num(cell2mat(p_test_values(:)))')>0);
    ind_train = (ismember(l,str2num(cell2mat(p_train_values(:)))')>0);
    
    train_GCF = GCF(ind_train);
    test_GCF = GCF(ind_test);
end

split_data_file = sprintf('split_GCF_%s',run_stamp);

s = whos('train_GCF');

if s.bytes > 2*1024^3
    disp 'the matrix train_GCF is huge'
    save(split_data_file,'train_GCF','test_GCF','ncents','-v7.3');%'cents',
else
    save(split_data_file,'train_GCF','test_GCF','ncents');%,'cents'
end


echo off
diary off
