function [GCF,dim_values] = LoadDataFiles(objects,angles,indeces,hights,scales,p_load_from_file)
% every file has a 2D matrix vid with complete sequence features. Columns
% store features of images. Matrix siz stores the sequence size
% (number of frames).

load('configuration');
%% preprocessing
% if  ~exist('p_load_from_file','var')
%     p_load_from_file=0;
% end
fileFolder = database.fileFolder;

dimension_names = {database.dimensions.name};
dim_values = cell(length(dimension_names));
i = strmatch('objects', dimension_names, 'exact');
if ~exist('objects','var') || strcmp(objects{1} , '*')    
    dim_values{i} = database.dimensions(i).values;
else
    dim_values{i} = objects;
end

i = strmatch('indeces', dimension_names, 'exact');
if ~exist('indeces','var') || strcmp(indeces{1} ,'*') 
    dim_values{i} = database.dimensions(i).values;
else
    dim_values{i} = indeces;
end

i = strmatch('angles', dimension_names, 'exact');
if ~exist('angles','var') || strcmp(angles{1}, '*')
    dim_values{i} = database.dimensions(i).values;
else
    dim_values{i} = strcmp;
end

i = strmatch('hights', dimension_names, 'exact');
if ~exist('hights','var') || strcmp(hights{1} ,'*')
    dim_values{i} = database.dimensions(i).values;
else
    dim_values{i} = hights;
end

i = strmatch('scales', dimension_names, 'exact');
if ~exist('scales','var') || strcmp(scales{1} ,'*')
    dim_values{i} = database.dimensions(i).values;
else
    dim_values{i} = scales;
end

%% if load from existing files

if exist('p_load_from_file','var') && p_load_from_file
    %tSeqs = [];
    GCF = query_GCF(dim_values);
    if ~isempty(GCF)
        return;
    end
end
%% save all data files first
%setupGaborKernel;%this is for Gabor filter

sz = [];
for i = 1: length(database.dimensions)
    sz(i) = length(database.dimensions(i).values);
end
if 1
tSizes = zeros(sz);
tImgs = cell(sz);

for i = 1:length(database.dimensions(1).values)
    folder = fullfile(fileFolder,database.dimensions(1).values{i});
    subfoldersdir = dir(folder);
    subfolders = {subfoldersdir.name};
    subfolders = subfolders(3:end);%to remove the . and .. folders from list
    for j = 1:length(subfolders)
        %m = regexp(subfolders{j},'.*_(\d*)','tokens');
        %index = str2double(m{1}{1});
        
        subfolder = fullfile(folder,subfolders{j});
        filesdir = dir(fullfile(subfolder,database.fileFilter));
        fileList = {filesdir.name};
        for file = fileList    
            m = regexp(file{1} , '.*\_A(\d)\_H(\d)\_S(\d)','tokens');
            x = str2double(m{1}{1});
            h= str2double(m{1}{2});
            s = str2double(m{1}{3});
            if ~validate_dim_val(h,database.dimensions(3).values)
                continue;
            end
            if ~validate_dim_val(s,database.dimensions(4).values)
                continue;
            end
            f = fullfile(subfolder,file{1});
            img = im2double(rgb2gray(imread(f)));
%             % masking the image
%             f_mask = fopen(fullfile(subfolder,'mask',strrep(file{1},'bmp','mask')));
%             v = fscanf(f_mask,'%d');
%             m = v(1); n = v(2);
%             mask = reshape(v(3:end),m,n);
%             mask(mask>0) = 1;
%             img = img.*mask;
            % extract features;
            tImgs{i,j,h,s,x} = featureExtraction(img,0);
            %tSizes(i,j,x,h,s) = size(tSeqs{i,k,s},2);%siz(3);
        end
    end
end
%clear fileList vid
filename = 'all_image_files_tensor';
save(filename, 'tImgs');
end
%% save data to file
if ~exist('tImgs','var')
    load('all_image_files_tensor');
end
tSeqs = cell(sz(logical(cell2mat({database.dimensions.category}))));
empty_seqs = [];
for i = 1:size(tImgs,1)
%     if ~validate_dim_val(i,database.dimensions(1).values)
%         continue;
%     end
    for j = 1:size(tImgs,2)
%         if ~validate_dim_val(j,database.dimensions(2).values)
%             continue;
%         end
        for h = 1:size(tImgs,3)
%             if ~validate_dim_val(h,database.dimensions(3).values)
%                 continue;
%             end
            for s = 1:size(tImgs,4)
%                 if ~validate_dim_val(s,database.dimensions(4).values)
%                     continue;
%                 end
                
                seq = squeeze(cell2mat(tImgs(i,j,h,s,:)));
                if ~isempty(seq)
                    tSeqs{i,j,h,s} = seq;
                else
                    empty_seqs = [empty_seqs;[i j h s]];
                end
                    
            end
        end
    end
end
clear tImgs;
% fill in holes in tSeqs
% for r = 1:size(empty_seqs,1)
%     if empty_seqs(r,3) == 2
%         pos = empty_seqs(r,:);
%         pos(3) = 3;
%         tSeqs{empty_seqs(r,:)} = tSeqs{pos};
%     else
%         if empty_seqs(r,3) == 3
%             pos = empty_seqs(r,:);
%             pos(3) = 2;
%             tSeqs{empty_seqs(r,:)} = tSeqs{pos};
%         end
%     end
% end

% save data
filename = 'all_data_files';
save(filename,'tSeqs');%,'tLetter','tSubject');
disp('learn mapping coeffiecient matrix for all datapoints');
GCF = getCoeffTensor_normalized(tSeqs,4,1);
%GCF = getCoeffTensor(tSeqs,1);
%end of preprocessing

%% return: take the part of interest
%tSeqs = tSeqs(tLetter,tRepeats,tSubject);
GCF = query_GCF(dim_values);
%Nseq=length(trainSeqs);%seq number for training

end