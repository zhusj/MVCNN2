mu1 = [2 3];
SIGMA1 = [1 1.5; 1.5 3];
r1 = mvnrnd(mu1,SIGMA1,100);
figure, plot(r1(:,1),r1(:,2),'+')
hold on;
mu2 = [5 4];
SIGMA2 = [2 1.5; 1.5 3];
r2 = mvnrnd(mu2,SIGMA2,100);
plot(r2(:,1),r2(:,2),'ro')


m1 = mean(r1);
m1 = repmat(m1,size(r1,1),1);
x1 = r1 - m1;
W1 = (x1'*x1)/100;

m2 = mean(r2);
m2 = repmat(m2,size(r2,1),1);
x2 = r2 - m1;
W2 = (x2'*x2)/100;

W = (W1+W2)/200;

r = [r1;r2];
m = mean(r);
m = repmat(m,size(r,1),1);
x = r - m;
T = (x'*x)/200;


P = inv(W)*T;

%[V,D]=eig(P)
[V,D]=eig(T,W)
%[U,S,V] = svd(P)

r_new = r*V;%(:,1);
figure,
%plot(r_new(1:100)',1,'b+');
plot(r_new(1:100,:),'b+');
hold on
%plot(r_new(101:200)',2,'ro');
plot(r_new(101:200,:),'ro');


