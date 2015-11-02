function get_average_test_result(file_expr,var_names,makeAvg)
if ~exist('makeAvg','var')
    makeAvg = [];
end

filelist = dir(file_expr);
filelist = {filelist.name};
if isempty(filelist)
    return;
end
total = cell(length(var_names),1);
scalarFlag = [];
for f = filelist
    s = load(f{1});
    for i = 1:length(var_names)
        l = s.(var_names{i});
        
        if isnumeric(l)==1
            scalarFlag(i) = true;
        else
            scalarFlag(i) = false;
        end
        
        %if ischar(l)
            
        if isstruct(l)
            
            if isempty(total{i})
                total{i} = l;
                continue;
            end
            vs = fieldnames(l);
            for j = 1:length(vs)
                if isnumeric(l.(vs{j}))
                    total{i}.(vs{j}) = [total{i}.(vs{j}),l.(vs{j})];
                end
            end 
        else
            l = l(:);
            total{i} = [total{i},l];
        end
    end

end

if isempty(makeAvg)
    makeAvg = mat2cell(scalarFlag(:),ones(length(scalarFlag),1));
end

avg = cell(length(total),1);
for i = 1:length(total)
    if ~makeAvg{i}
        avg{i} = total{i};
    else
        
        if isstruct(total{i})
            vs = fieldnames(total{i});
            for j = 1:length(vs)
                if isnumeric(total{i}.(vs{j}))
                    avg{i}.(vs{j}) = sum(total{i}.(vs{j}),2)/size(total{i}.(vs{j}),2);
                end
            end
        else
            avg{i} = sum(total{i},2)/size(total{i},2);
        end
    end
end
s_out = cell2struct(avg,var_names,1);

new_name = strrep(file_expr,'*','overall');
save(new_name,'-struct','s_out');