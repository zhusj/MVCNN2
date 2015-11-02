seq = tSeqs{26,1,5};

ncent = size(seq,2);% mapping centers, Need to know if this number affects the performance???????
icent=((1:ncent)-1)'*2*pi/ncent;
cent = [cos(icent) sin(icent)];
g =full(directs_grbf_reg(cent,cent));


A = seq';
K = A*A';
[Y,B,T,U,W,C,S]=KPLS(A,K,g,3);

figure,plot3(T(:,1),T(:,2),T(:,3));
hold on
for i=1:3:size(A,1)
    v = A(i,:);
    v1 = v*W;
    %v1 = v1./S;
    plot3(v1(1),v1(2),v1(3),'r*');
    v2 = v*B.M;
    plot3(v2(1),v2(2),v2(3),'b*');
end