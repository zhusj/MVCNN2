f ='c:\Amr\ChanSu\MyCode\avletters\Lips\Z1_Anya-lips.mat';
load(f);
i = strfind(f,'.mat');
i = max(i);
p = f(1:i-1);
mkdir(p);
for i = 1:siz(3)
    c = reshape(vid(:,i),siz(1),siz(2));
    colormap(gray);
    c = c/max(max(c));
    imwrite(c,sprintf('%s/image%d.jpg',p,i));
end