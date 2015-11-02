function delete_files_expr(file_expr,varargin)
filelist = dir(file_expr);
filelist = {filelist.name};
if isempty(filelist)
    return;
end

for f = filelist
    tic
    delete(f{1});
    toc
end

% for f = filelist
%     s = rmfield(load(f{1}),varargin(:));
%     save(f{1},'-struct','s');
% end
