function remove_field_files(file_expr,varargin)
filelist = dir(file_expr);
filelist = {filelist.name};
if isempty(filelist)
    return;
end

% if nargin<2% ~exist(var_names,'var')
%     s = load(filelist{1});
%     var_names = fieldnames(s);
% end
A_train = [];
B_train = [];
for f = filelist
    f{1}
    tic
    %save(f{1},'A_train','-append');
    save(f{1},'A_train','B_train','-append');
    toc
end

% for f = filelist
%     s = rmfield(load(f{1}),varargin(:));
%     save(f{1},'-struct','s');
% end
