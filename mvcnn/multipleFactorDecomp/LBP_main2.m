function I2=LBP_main2(I,Nc,Nr) 
global lbp_mapping;
% SP=[-1 -1; -1 0; -1 1; 0 -1; -0 1; 1 -1; 1 0; 1 1];
% I2=lbp(I,SP,0,'i'); %LBP code image using sampling points in SP
%                   %and no mapping. Now H2 is equal to histogram
%                   %of I2.
%figure;
%colormap(gray);
if ~exist('lbp_mapping','var')
    lbp_mapping=getmapping(8,'u2'); 
    %load('configuration',)
end
I2 = [];
%H = fspecial('gaussian',[3 3],0.5);

%Nc = 4;
%Nr = 3;
[r,c] = size(I);
cols = floor((0:Nc)/Nc*c);
rows = floor((0:Nr)/Nr*r);
% for i = 1:N
%     I1 = I(:,cols(i)+1:cols(i+1));
%     %t = lbp(I,3,8,lbp_mapping,'nh');
%     t = lbp(I1,3,8,lbp_mapping,'h');
%     I2 = [I2;t(:)];
% end
% return;
for j = 1:Nr
for i = 1:Nc
    I1 = I(rows(j)+1:rows(j+1),cols(i)+1:cols(i+1));

    for n = 1:8%[0,1,2,3]
        %subplot(2,2,n);
        %t = lbp(I,2^n,8,mapping,'nh');
        t = lbp(I1,n,8,lbp_mapping,'nh');
        
        t = t(2:end-1);%remove the extremes, usually extreme bins includes unwanted informatoin
        %t = imfilter(t,H,'replicate');
        %imagesc(t);
        I2 = [I2;t(:)];
        %I2 = t;
    end
end
end
%I2 = I2(:);