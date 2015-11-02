function F = featureExtraction2(A,sequence,m,n)%this function is devoted to return the features corresponding to sequence/single image
%^sequence is boolean for identifying. m,n are size of frame if sequence is
if nargin<2
    sequence = 1;
    m=80;n=60;
end

%true
%% As Is
 %F = A;
len = size(A,2);
pts = ((1:len)-(len/2.0))/(len/2.0);
l = normpdf(pts,0,0.5);
l = repmat(l,size(A,1),1);
F = A.*l;
%  perc = floor(size(A,2)/7.0);
%  F = A(:,perc:end-perc);
  m = max(max(F));
  F = F/m;
 return;
 %  m = max(F);
%  m = repmat(m,size(F,1),1);
%  F = F./m;
% % [U,S,V]= svd(A,0);
% % F = V;

%% 2D - DCT
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

% %% Gabor filter
% numframes = min(size(A));
% C = reshape(A,m,n,numframes);
% d = Gabor(C(:,:,1));
% F = zeros(length(d),numframes);
% F(:,1) = d;
% if sequence
%     for i=2:numframes
%         F(:,i) = Gabor(C(:,:,i));
%     end
% else
%     F = dct2(A);
% end
% return;
%% 2D - HoG
numframes = min(size(A));
C = reshape(A,m,n,numframes);
d = HoG(C(:,:,1));
F = zeros(length(d),numframes);
F(:,1) = d;
if sequence
    for i=2:numframes
        F(:,i) = HoG(C(:,:,i));
    end
else
    F = dct2(A);
end
return;
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