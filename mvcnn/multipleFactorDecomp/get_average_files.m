function get_average_files(file_expr,var_names)
filelist = dir(file_expr);
filelist = {filelist.name};
if isempty(filelist)
    return;
end


if nargin<2% ~exist(var_names,'var')
    var_names = {'recog_rate','recog_rate_unit', 'unit_conflict'};
else
    if strcmp(var_names,'*')
        s = load(filelist{1});
        var_names = fieldnames(s);
    end
end

var_names_vrnc = var_names;
for i = 1:length(var_names)
    var_names_vrnc{i}  = [var_names{i},'_variance'];
end

total = cell(length(var_names),1);
variance = cell(length(var_names),1);
for f = filelist
    s = load(f{1});
    for i = 1:length(var_names)
        l = s.(var_names{i});
        if ~isnumeric(l)
            continue;
        end
        if isempty(total{i})
            
            total{i} = l;
            variance{i} = l.^2;
        else
            if sum(size(l) ~= size(total{i}))
                continue;
            end
            total{i} = total{i}+l;
            variance{i} = variance{i}+l.^2;
        end
    end
end

n = length(filelist);

avg = cell(length(total),1);
for i = 1:length(total)
    avg{i} = total{i}/n;
    variance{i} = (variance{i} -  n*(avg{i}).^2)/(n-1);%unbiased estimate of sample variance
end

new_name = strrep(file_expr,'*','overall');

s_out = cell2struct(avg,var_names,1);
save(new_name,'-struct','s_out');

s_out = cell2struct(variance,var_names_vrnc,1);
save(new_name,'-struct','s_out','-append');