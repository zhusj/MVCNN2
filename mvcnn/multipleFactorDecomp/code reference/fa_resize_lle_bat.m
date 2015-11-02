
subject_name = 's4_bigFakesmile_04_bmp';
% subject_name = 's1_SoftAffection_02_bmp';
pathname = 'H:\Project\multifaceanimation\data\texture\';
outfile_name = strcat(pathname,subject_name,'-normalized')
load(outfile_name);%,'croppedImage','croppedParam');    


FrameNUM = size(croppedImage,2);
RESIZE1 = 48;
RESIZE2 = 96;
resizedImage48 = zeros(FrameNUM,RESIZE1*RESIZE1);
resizedImage96 = zeros(FrameNUM, RESIZE2*RESIZE2);

for i=1:FrameNUM

    Img_i = croppedImage{i};
    resizedImage48(i,:) = reshape(imresize(Img_i,[RESIZE1 RESIZE1]),1,RESIZE1*RESIZE1);
    resizedImage96(i,:) = reshape(imresize(Img_i,[RESIZE2 RESIZE2]),1, RESIZE2*RESIZE2);
end

newfilename1 = strcat(subject_name,'resized48');
save(newfilename1,'resizedImage48');

Y = lle(resizedImage48(1:480,:)',55,4)
figure
plot3(Y(1,:), Y(2,:), Y(3,:),'r.-')

show_batches1(resizedImage48(40:40:480,:),48, 48, 12,0)


figure
for i=40:40:480
    frame_str = int2str(i);
    text(Y(1,i), Y(2,i), Y(3,i),frame_str)
end
hold on
plot3(Y(1,:), Y(2,:), Y(3,:),'r.-')

% % newfilename1 = strcat(subject_name,'resized96');
% % save(newfilename1,'resizedImage96');
% 
% for k=15:5:80
%     Y = lle(resizedImage48',k,4)
% 	figure
% 	plot3(Y(1,:), Y(2,:), Y(3,:),'r.-')
%     title(['K=' int2str(k)])
%     drawnow
% end
% 
% %   
% % % for k=80:20:300
% % %     Y = lle(resizedImage48',k,4)
% % % 	figure
% % % 	plot3(Y(1,:), Y(2,:), Y(3,:),'r.-')
% % % %     plot(Y(1,:), Y(2,:),'r.-')
% % %     title(['K=' int2str(k)])
% % %     drawnow
% % % end
% % 
