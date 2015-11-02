function I2=LBP_main(I) 
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
for n = 1:8%[0,1,2,3]
    %subplot(2,2,n);
    %t = lbp(I,2^n,8,mapping,'nh');
    t = lbp(I,n,8,lbp_mapping,'nh');
    
    %t = imfilter(t,H,'replicate');
    %imagesc(t);
    I2 = [I2;t(:)];
    %I2 = t;
end
%I2 = I2(:);