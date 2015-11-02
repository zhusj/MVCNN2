function [F_,oP1_] = featureExtraction(featureType,A,sequence,Nr,Nc,m,n)%this function is devoted to return the features corresponding to sequence/single image
%^sequence is boolean for identifying. m,n are size of frame if sequence is
if ~exist('sequence','var')
    sequence = 1;
end
if ~exist('m','var')
    m=60;n=80;
end
if ~exist('Nc','var')
    Nc = 1; Nr = 1;
end

%true

if ~exist('sequence','var') || sequence
    %F = A;
     perc = floor(size(A,2)/7.0);
     perc = 0;
     F_ = A(:,1:end-perc);
     mx = double(max(max(F_)));
     %F = F/mx;
    numframes = min(size(F_)); 
    C = reshape(F_,m,n,[]);
else
    [m n l] = size(A);
    numframes = 1;
    if l>1 && ~(strcmp(featureType,'RGB')||strcmp(featureType,'rgbHOG')) 
        %the image is RGB then convert it
        C = rgb2gray(A);
    else
        C = A;
    end
    %C=A;
end

 
%  m = max(F);
%  m = repmat(m,size(F,1),1);
%  F = F./m;
% % [U,S,V]= svd(A,0);
% % F = V;
% return;

%% subsample
if(strcmp(featureType,'RGB'))
    R = C(:,:,1);
    G = C(:,:,2);
    B = C(:,:,3);
    
    %[r,c] = size(I);
    cols = floor((0:Nc)/Nc*n);  
    rows = floor((0:Nr)/Nr*m);

%     r_step = floor(m/Nr);
%     c_step = floor(n/Nc);
    r_hist = [];g_hist = [];b_hist = [];
    x = 0:0.1:1;
    
    for j = 1:Nr
    for i = 1:Nc
            h = hist(reshape(R(rows(j)+1:rows(j+1),cols(i)+1:cols(i+1)),1,[]),x);
            r_hist = [r_hist,h(2:end)];
            h = hist(reshape(G(rows(j)+1:rows(j+1),cols(i)+1:cols(i+1)),1,[]),x);
            g_hist = [g_hist,h(2:end)];
            h = hist(reshape(B(rows(j)+1:rows(j+1),cols(i)+1:cols(i+1)),1,[]),x);
            b_hist = [b_hist,h(2:end)];
    end
    end
    r_hist = r_hist/sum(r_hist);
    g_hist = g_hist/sum(g_hist);
    b_hist = b_hist/sum(b_hist);
    %K = filter2(fspecial('average',3),A)/255;
    F_ = [r_hist',g_hist',b_hist'];
%L = medfilt2(J,[3 3]);
%conv2
end

%% subsample
if(strcmp(featureType,'RESIZE'))
    %K = filter2(fspecial('average',3),A)/255;
    F_ = imresize(C,[Nr Nc]);
%L = medfilt2(J,[3 3]);
%conv2
end
%% LBP
if(strcmp(featureType,'LBP'))
    d = LBP_main2(C(:,:,1),Nr,Nc);
    F_ = zeros(length(d),numframes);
    F_(:,1) = d;
    for i=2:numframes
        F_(:,i) = LBP_main2(C(:,:,i),Nr,Nc);
    end
    return;
end



%% 2D - DCT
if(strcmp(featureType,'DCT'))
    F_ = zeros(size(A));
    if sequence
        for i=1:size(A,2)
            c = reshape(A(:,i),m,n);
            d = dct2(c);
            F_(:,i) = reshape(d,numel(d),1);
        end
    else
        F_ = dct2(A);
    end
    return;
end
%% Gabor filter
if(strcmp(featureType,'GABOR'))
    d = Gabor(C(:,:,1));
    F_ = zeros(length(d),numframes);
    F_(:,1) = d;
    for i=2:numframes
        F_(:,i) = Gabor(C(:,:,i));
    end
    return;
end

%% 2D - HoG
if(strcmp(featureType,'HOG'))
    %adjust image size so that min cell size is hold and preserve aspect ratio
    min_cell_sz = 3; max_cell_sz = 15;
    min_ratio = min([m n]./[Nr Nc]);
        
    F_ = [];
    for winsz_ratio = 1%[4,10,20]
        img = C(:,:,1);
        if min_ratio < min_cell_sz
            img = imresize(img,min_cell_sz/min_ratio);
        else
            if min_ratio > max_cell_sz
                img = imresize(img,max_cell_sz/min_ratio);
            end
        end
        
        [d,oP1_]= HOG(img,Nr,Nc,16);
        
        F = zeros(length(d),numframes);
        F(:,1) = d;

        for i=2:numframes
            F(:,i) = HOG(C(:,:,i),winsz_ratio);
        end
        F_ = [F_;F];
    end
    return;
end
%% 2D - hpHOG
if(strcmp(featureType,'rgbHOG'))
    %adjust image size so that min cell size is hold and preserve aspect ratio
    min_cell_sz = 4; max_cell_sz = 10;
    min_ratio = min([m n]./[Nr Nc]);
        
    F_ = [];
    for winsz_ratio = 1%[4,10,20]
        %img = C(:,:,1);
%         if min_ratio < min_cell_sz
%             C = imresize(C,min_cell_sz/min_ratio);
%         else
%             if min_ratio > max_cell_sz
%                 C = imresize(C,max_cell_sz/min_ratio);
%             end
%         end
        
        %[d,oP1_]= HOG(img,Nr,Nc,9);
        d= HOG(C(:,:,1),Nr,Nc,9);
        %F = d;
        F = zeros(length(d),numframes);
        F(:,1) = d;

        for i=2:size(C,3)
            F(:,i) = HOG(C(:,:,i),Nr,Nc,9);
        end
        F_ = [F_;F];
    end
    return;
end
%% 2D - SIFT
%     F = zeros(size(A));
%     if sequence
%         for i=1:size(A,2)
%             c = reshape(A(:,i),m,n);
%             d = dct2(c);
%             F(:,i) = reshape(d,numel(d),1);
%         end
%     else
%         F = dct2(A);
%     end
% %%2D - Zernike moments
%     F = zeros(size(A));
%     if sequence
%         for i=1:size(A,2)
%             c = reshape(A(:,i),m,n);
%             d = dct2(c);
%             F(:,i) = reshape(d,numel(d),1);
%         end
%     else
%         F = dct2(A);
%     end
end