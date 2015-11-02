function str = cell2string(c)
    str = [];
    for i = 1:length(c)
        str = [str,'-',c{i}];
    end
    
    str = str(2:end);
end