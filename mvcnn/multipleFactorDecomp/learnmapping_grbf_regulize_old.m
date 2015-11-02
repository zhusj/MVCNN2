% learn GRBF mapping
% input:
% D : original data matrix Nxd  N: is the number of samples, d: is the dimensionality of the input images
% X : embedding space representation. Nxe matrix N: number of samples, e: the dimesnionality of the embedding space
% cent: centers: exm e: the dimesnionality of the embedding space

function [CF,cent]=learnmapping_grbf_regulize(p_Data,cent,points)
%CF = learnmapping_KPLS(D);
%return

%if centers number is not given then set it (better)
if ~exist('cent','var') || isempty(cent)
    ncent = 7;% mapping centers, Need to know if this number affects the performance???????
    icent=((1:ncent)-1)'*2*pi/ncent;
    cent = [cos(icent) sin(icent)];
    %Nb=length(cent);
    %d=2;
end
lambda = -0.1;
thresh = 1E-5;

D = p_Data;
%%%%%
%cent = cent(2:end,:);
%%%%%
count = 0;
% while true

if ~exist('points','var')||isempty(points)
    % embed on a circle
    % % %measure how close every two consecutive frames:
    % % d = dist2(D,D);
    % % d = [diag(d,1);d(1,end)];
    % indx = 1:size(D,1);
    % indx2 = [indx(2:end),1];
    % d = row_wise_dist(D(indx,:) ,D(indx2,:),1,1);
    % d = d/sum(d);
    % %accumlate distances so that last value is 1:
    % accumlate = tril(ones(length(d),length(d)),0);
    % d = accumlate*d;
    % d = [0;d(1:end-1)];
    %%%
    frames=size(D,1);%angletrain{seq};
    d = ((1:frames)'-1)/frames;
    %get embedding on circle.
    %d = normrnd(d,.005);
    t=d*2*pi;
else
    t = points(:);
end
[t,ord] = sort(t);
D = p_Data(ord,:);
X=[cos(t) sin(t)];
%implementation of higher dimensional regularization
%cf = (G'G+lambdag)-1 G'y

N=size(D,2);
G =directs_grbf_reg(X,cent);
g =directs_grbf_reg(cent,cent);
% G =directs_grbf_reg(t,icent);
% g =directs_grbf_reg(icent,icent);

F =G'*D;
B =full(G'*G-lambda*g);
CF = B\F;

count = count+1;

% if exist('CF_old','var')
%     diff = norm(CF-CF_old);
%     if diff<thresh || count>20
%         break;
%     end
% end
%CF_old = CF;
% D = G*CF;
% disp 'The reconstruction error is:'
% err = norm(D-p_Data)
%err = mean(row_wise_dist(D',p_Data','corr',1,1))
% if err <= 0.01
%     break;
% end

%lambda = 0;
%lambda = lambda/2;
%%%%%
% CF = CF(2:end,:);
%CF = inv(G'*G-lambda*g)*F;
end
%count
