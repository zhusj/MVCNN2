function [GCF,tSeqs,tLetter,tSubject] = LoadDataAVLetterFiles(letters,subjects,repeats,p_load_from_file)
% every file has a 2D matrix vid with complete sequence features. Columns
% store features of images. Matrix siz stores the sequence size
% (number of frames).
fileFolder = fullfile('C:','Amr','ChanSu','MyCode','avletters','Lips');
%fileFolder = fullfile('C:','Amr','ChanSu','MyCode','avletters','Gabor','setting2');
tSizes = [];
tLetter = [];
%% aaa
% if  ~exist('filename','var')
%     filename='test_data';
% end

if strcmp(letters{1} , '*') ||  ~exist('letters','var')
    letters = {'A';'B';'C';'D';'E';'F';'G';'H';'I';'J';'K';'L';'M';'N';'O';'P';'Q';'R';'S';'T';'U';'V';'W';'X';'Y';'Z'};
end
if strcmp(subjects{1}, '*') ||  ~exist('subjects','var')
    subjects = {'Anya';'Bill';'Faye';'John';'Kate';'Nicola';'Stephen';'Steve';'Verity';'Yi'};
end
if strcmp(repeats{1} ,'*') ||  ~exist('repeats','var')
    repeats = {'1';'2';'3'};
end
%%

tLetter = zeros(length(letters),1);
for i = 1:length(letters)
        tLetter(i) = (letters{i}-'A'+1);
end

tRepeats = zeros(length(repeats),1);
for i = 1:length(repeats)
        tRepeats(i) = str2double(repeats{i});
end
tSubject = zeros(length(subjects),1);
for i = 1:length(subjects)
    tSubject(i) = getSubjectLabel(subjects{i});
end

%filename = 'all_data_files';
if exist('p_load_from_file','var') && p_load_from_file
    load('all_GCF');
    GCF = GCF(tLetter,tRepeats,tSubject,:,:);
    return;
end

filename = 'all_data_files';
save(filename,'tSeqs');%,'tLetter','tSubject');
disp('learn mapping coeffiecient matrix for all datapoints');
GCF = getCoeffTensor(tSeqs);
filename = 'all_GCF';
save(filename,'GCF');
%end of preprocessing

%take the part of interest
GCF = GCF(tLetter,tRepeats,tSubject,:);
%Nseq=length(trainSeqs);%seq number for training

end