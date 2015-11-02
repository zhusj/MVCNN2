function T_part = query_tensor(T,p_names,p_values)
%dim_values is a cell array, one cell for every dimension sorted by
%dimension order

if mod(nargin,2)<1 %p_namex and p_valuesx should comes in pairs
    T_part = [];
    return;
end

if length(p_names) ~= length(p_values)
   T_part = [];
    return;
end 

%if ~exist('all_GCF.mat','file') || ~exist('configuration.mat','file')
if ~exist('configuration.mat','file')
    T_part = [];
    return;
end
load('configuration','database');
%load('all_GCF');
%dim_num = sum(cell2mat({database.dimensions.category}));

dimension_names = {database.dimensions.name};

for i = 1:length(p_names)
    indx(i) = strmatch(p_names{i}, dimension_names, 'exact');
end

dim_num = length(size(T));
S = 'T(';
v = [];
for i = 1:dim_num
    c = find(indx==i);
    if isempty(c)
        S = [S,':',','];
    else
        v = getClassLabel(p_values{c},database.dimensions(i).values);
        S = [S,mat2str(v'),','];
    end
end
S = [S(1:end-1),')'];
T_part = eval(S);

% switch dim_num
%     case 2
%         GCF_part = GCF(V{1},V{2},:);
%     case 3
%         GCF_part = GCF(V{1},V{2},V{3},:);
%     case 4
%         GCF_part = GCF(V{1},V{2},V{3},V{4},:);
%     case 5
%         GCF_part = GCF(V{1},V{2},V{3},V{4},V{5},:);
%     case 6
%         GCF_part = GCF(V{1},V{2},V{3},V{4},V{5},V{6},:);
%     case 7
%         GCF_part = GCF(V{1},V{2},V{3},V{4},V{5},V{6},V{7},:);
% end