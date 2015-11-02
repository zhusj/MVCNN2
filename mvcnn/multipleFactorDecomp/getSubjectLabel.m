function l = getSubjectLabel(s)
    %1  Anya
    %2  Bill
    %3  Faye
    %4  John
    %5  Kate
    %6  Nicola
    %7  Stephen
    %8  Steve
    %9  Verity
    %10 Yi
  subjects = {'Anya','Bill','Faye','John','Kate','Nicola','Stephen','Steve','Verity','Yi'};
  %f1 = strfind(s,'_');
  %f2 = strfind(s,'-');
  %str = s(f1+1:f2-1);
  %l = strmatch(str, subjects, 'exact');
  if s=='*'
    l=[1:length(subjects)]';
  else
    l = strmatch(s, subjects, 'exact');
  end
  
    