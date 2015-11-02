function l = getClassLabel(s,list)
  if strcmp(s,'*')
    l=(1:length(list))';
  else
    %l = strmatch(s, list, 'exact');
    for i = 1:length(s)
        l(i) = strmatch(s{i}, list, 'exact');
    end
  end
  
  if isempty(l)
      l = -1;
  else
      l = l(:);
  end
    