function l = getLetterName(idx,b)
if nargin < 2
    b = 26;
end

%% Groups
%      'A','K','L';
%      'B','P';
%      'C','D','E','T','Z';
%      'F';
%      'G','J';
%      'H';
%      'I','R';
%      'M';
%      'N';
%      'O','Q','U';
%      'S','X';
%      'V';
%      'W';
%      'Y'

%     LIST = [
%         'A'
%         'B'
%         'C'
%         'C'%'D'
%         'C'%'E'
%         'F'
%         'G'
%         'H'
%         'I'
%         'G'%'J'
%         'A'%'K'
%         'A'%'L'
%         'M'
%         'N'
%         'O'
%         'B'%'P'
%         'O'%'Q'
%         'I'%'R'
%         'S'
%         'C'%'T'
%         'O'%'U'
%         'V'
%         'W'
%         'S'%'X'
%         'Y'
%         'C'%'Z'
%  ];
%%
LIST = [
        'A'
        'B'
        'C'
        'D'
        'E'
        'F'
        'G'
        'H'
        'I'
        'J'
        'K'
        'L'
        'M'
        'N'
        'O'
        'P'
        'Q'
        'R'
        'S'
        'T'
        'U'
        'V'
        'W'
        'X'
        'Y'
        'Z'
        ];
    %l = char('A'+ceil(double(idx)/reptn)-1);
    l = LIST(mod(idx-1,b)+1);