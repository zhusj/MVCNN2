function l = getSubjectName(s)
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
  if s == '*'
      l=s;
      return;
  end
  subjects = {'Anya','Bill','Faye','John','Kate','Nicola','Stephen','Steve','Verity','Yi'};
  
  if ischar(s)
      s = str2double(s);
  end
  
%   if s<1 || s>10
%       l = 'None';
%       return;
%   end
  l = subjects{s};
  
    