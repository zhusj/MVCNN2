function l = validate_dim_val(s,list)
  if isempty(list)
      l = 1;
      return;
  end
  
  if isnumeric(s)
      s = sprintf('%d',s);
  end
  l = strmatch(s, list, 'exact');
  if isempty(l)
      l = 0;
  end
    