function I2=LBP_main_mask(I,mask) 
global lbp_mapping;
% SP=[-1 -1; -1 0; -1 1; 0 -1; -0 1; 1 -1; 1 0; 1 1];
% I2=lbp(I,SP,0,'i'); %LBP code image using sampling points in SP
%                   %and no mapping. Now H2 is equal to histogram
%                   %of I2.
%figure;
%colormap(gray);
if ~exist('lbp_mapping','var')
    lbp_mapping=getmapping(8,'u2'); 
end
I2 = [];
%H = fspecial('gaussian',[3 3],0.5);
figure;
for n = 1:8%[0,1,2,3]
    %t = lbp(I,n,8,lbp_mapping,'nh');
    t = lbp(I,n,8,lbp_mapping,'i');
    diff = size(mask)-size(t);%get the difference in size between mask and result of lbp
    mask = mask(diff(1)/2:end-diff(1)/2-1,diff(2)/2:end-diff(2)/2-1);
    subplot(1,2,1), imagesc(t);
    t = t.*mask;
    subplot(1,2,2),imagesc(t);
    t=hist(t(:),0:(lbp_mapping.num-1));
    t = t(2:end);
    t=t/sum(t);
    
    %t = imfilter(t,H,'replicate');
    %imagesc(t);
    I2 = [I2;t(:)];
    %I2 = t;
end
%I2 = I2(:);