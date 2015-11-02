
subject_name = 's4_bigFakesmile_04_bmp';
pathname = 'H:\Project\multifaceanimation\data\texture\';
outfile_name = strcat(pathname,subject_name,'-normalized')
load(outfile_name);%,'croppedImage','croppedParam');    

filename='list.txt';
[fname]=textread(strcat(pathname,subject_name,'\',filename),'%s');

        





FrameNUM = size(fname,1);

IWIDTH = 48;
IHEIGHT = 48;
DATA_DIM = IWIDTH*IHEIGHT;


% HISTOGRAM_THRESHOLD_NUM = 15;
HISTOGRAM_THRESHOLD_NUM = 8;
NORMAL_THRESHOLD_VALUE = 1.00;

histogram_scaled_threshold_image = zeros(FrameNUM, DATA_DIM);

resizedImage = zeros(FrameNUM, DATA_DIM);
for k=1:FrameNUM
    I2_rgb =imread(strcat(pathname,subject_name,'\',fname{k}));
    Img_k = rgb2gray(I2_rgb);
    resizedImage(k,:) = double(reshape(imresize(Img_k, [IWIDTH IHEIGHT]),1, DATA_DIM));
%     histogram_dist = hist(resizedImage(k,:)',256);
%     histogram_scaled_image = resizedImage(k,:)/max(find(histogram_dist>=HISTOGRAM_THRESHOLD_NUM));
%     threshold_image = histogram_scaled_image;
%     threshold_image(histogram_scaled_image>1) = 1;
%     histogram_scaled_threshold_image(k,:)=threshold_image;
end

 Y = lle(resizedImage(1:500,:)',55,4)
figure
plot3(Y(1,:), Y(2,:), Y(3,:),'r.-')


return

for k=5:5:100
    Y = lle(resizedImage(1:500,:)',k,4)
	figure
	plot3(Y(1,:), Y(2,:), Y(3,:),'r.-')
%     plot(Y(1,:), Y(2,:),'r.-')
    title(['K=' int2str(k)])
    drawnow
end

