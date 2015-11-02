mu1 = [2 3];
SIGMA1 = [1 1.5; 1.5 3];
r1 = mvnrnd(mu1,SIGMA1,100);
figure, subplot(2,2,1);
plot(r1(:,1),r1(:,2),'+')
hold on;
mu2 = [5 0];
SIGMA2 = [3 0; 0 4];
r2 = mvnrnd(mu2,SIGMA2,100);
plot(r2(:,1),r2(:,2),'ro')

mu3 = [4 1];
SIGMA3 = [3 1.5; 1.5 1];
r3 = mvnrnd(mu3,SIGMA3,100);
plot(r3(:,1),r3(:,2),'cx')
r = [r1;r2;r3];

y = [ones(length(r1),1);2*ones(length(r2),1);3*ones(length(r3),1)];

mr =mean(r);
r = bsxfun(@minus,r,mr);

my = mean(y);
y = y-my;

%% PLS
w = r'*y;
w = w / (w'*w);
%hold, line([0,w(1)],[0,w(2)]);
r_new = r*(w*w');
%% Kernel PLS
I = eye(length(unique(y)));
lbl = [];
for i = unique(y')
    indx = find(y==i);
    lbl =  [lbl;repmat(I(i+2,:),length(indx),1)];
end
n = size(r,1);
A=exp(-1*pdist2(r,r));%full(directs_grbf_reg(Y,Y));
% r=dist2(X,X);
% A = exp(-0.5*r);
% clear r;
K = A*A';%get_gram_kernel(A,A);
w = ones(n,1);
w= w*w';
w = eye(n)-w/n;
K = w*K*w';%normalizing K using equation (13) in the paper
[lbl_new,B,T,U] = KPLS(A,K,lbl);
%[~,r_new] = max(lbl_new,[],2);
r_new = lbl_new;
% 
% figure,
% title('T');
% plot3(T(1:100,1),T(1:100,2),T(1:100,3),'b+');
% hold on
% plot3(T(101:200,1),T(101:200,2),T(101:200,3),'ro');
% plot3(T(201:300,1),T(201:300,2),T(201:300,3),'cx');
% 
figure,
title('U');
plot3(U(1:100,1),U(1:100,2),U(1:100,3),'b+');
hold on
plot3(U(101:200,1),U(101:200,2),U(101:200,3),'ro');
plot3(U(201:300,1),U(201:300,2),U(201:300,3),'cx');

W = ((K*U));
figure,
title('W');
% plot3(W(:,1),W(:,2),W(:,3),'r.');
% W = W';
% figure,
% title('W');
% plot3(W(:,1),W(:,2),W(:,3),'b.');

plot3(W(1:100,1),W(1:100,2),W(1:100,3),'b+');
hold on
plot3(W(101:200,1),W(101:200,2),W(101:200,3),'ro');
plot3(W(201:300,1),W(201:300,2),W(201:300,3),'cx');

%% plotting
 subplot(2,2,2)
% plot(r_new(1:100,1),r_new(1:100,2),'b+');
% hold on
% plot(r_new(101:200,1),r_new(101:200,2),'ro');
% plot(r_new(201:300,1),r_new(201:300,2),'cx');

subplot(2,2,3)
plot(r_new(1:100,1),1,'b+');
hold on
plot(r_new(101:200,1),2,'ro');
plot(r_new(201:300,1),3,'cx');

subplot(2,2,4)
plot(r_new(1:100,2)',1,'b+');
hold on
plot(r_new(101:200,2)',2,'ro');
plot(r_new(201:300,2),3,'cx');
return
%%
 GCF = K;
DDT = GCF*GCF';
%         [U, SS, UT] = svd(DDT);
    [U, SS, UT] = svds(DDT, size(GCF,1));
    S = sqrt(SS);
    invS = inv(S);
    U = invS*U'*GCF;
plotting(U',y,0);
