function T_part = query_tensor(T,p_name1,p_values1,p_name2,p_values2,p_name3,p_values3)
%dim_values is a cell array, one cell for every dimension sorted by
%dimension order

if mod(nargin,2)<1 %p_namex and p_valuesx should comes in pairs
    T_part = [];
    return;
end

%if ~exist('all_GCF.mat','file') || ~exist('configuration.mat','file')
if ~exist('configuration.mat','file')
    T_part = [];
    return;
end
load('configuration',database);
%load('all_GCF');
%dim_num = sum(cell2mat({database.dimensions.category}));

dimension_names = {database.dimensions.name};

indx = strmatch(dimension_name, dimension_names, 'exact');

dim_num = length(size(T));
S = 'T(';
v = [];
for i = 1:dim_num 
    v = getClassLabel(dim_values{i},database.dimensions(i).values);
    l = length(unique(v));
    if  (l == length(v)) && (l == length(database.dimensions(i).values))
        S = [S,':',','];
    else
        S = [S,mat2str(v'),','];
    end
end
S = [S,':)'];
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