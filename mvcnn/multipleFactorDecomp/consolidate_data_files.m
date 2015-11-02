function file_new_name = consolidate_data_files(file_expr,var_names)
filelist = dir(file_expr);
filelist = {filelist.name};
if isempty(filelist)
    return;
end

if ~exist('var_names','var')||isempty(var_names)
    s = load(filelist{1});
    var_names = fieldnames(s);
end

total = cell(length(var_names),1);
for f = filelist
    s = load(f{1});
    for i = 1:length(var_names)
        l = s.(var_names{i});
        %l = l(:)';
        total{i} = [total{i};l];
    end

end


% avg = cell(length(total),1);
% for i = 1:length(total)
%     avg{i} = sum(total{i})/size(total{i},1);
% end
s_out = cell2struct(total,var_names,1);


file_new_name = strrep(file_expr,'*','ALL');
save(file_new_name,'-struct','s_out');