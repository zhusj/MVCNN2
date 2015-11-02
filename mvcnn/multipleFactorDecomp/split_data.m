function [split_data_file] = split_data(file_name,var_name,dimension_name,p_values,run_stamp)
%clear
%%preprocessing
echo on
S = load(file_name,var_name);
T = getfield(S,var_name);
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
if iscell(T) && (numel(size(T))>2)
    new_T = query_tensor(T,{dimension_name},{p_values});
    %test_GCF = query_tensor(T,{dimension_name},{p_test_values});
end
if isstruct(T)
    l = cell2mat({T.(dimension_name)});
    ind = (ismember(l,str2num(cell2mat(p_values(:)))')>0);
    new_T = T(ind);
end

if exist('run_stamp','var')
    split_data_file = sprintf('split_%s_%s',var_name,run_stamp);
else
    split_data_file = sprintf('split_%s_%s',var_name,dimension_name);
end


S_out = struct;
%setfield(S_out,var_name,new_T);
S_out = cell2struct({new_T} ,{var_name},1);



    

s = whos('new_T');

if s.bytes > 2*1024^3
    disp 'the matrix train_GCF is huge'
    save(split_data_file ,'-struct','S_out','-v7.3');
else
    save(split_data_file ,'-struct','S_out');
end

