% learn GRBF mapping
% input:
% D : original data matrix Nxd  N: is the number of samples, d: is the dimensionality of the input images
% X : embedding space representation. Nxe matrix N: number of samples, e: the dimesnionality of the embedding space
% cent: centers: exm e: the dimesnionality of the embedding space

function CF=learnmapping_latent_CF(p_input,p_output,p_labels)
%I have input points, output points and need regression matrix.

%compute centers as the average of latent points for every class
cent = zeros(numel(unique(p_labels)),size(p_input,2));
c = 1;
for i = p_labels'
    indx = (p_labels==i);
    l = p_input(indx,:);
    cent(c,:) = mean(l);
    c = c+1;
end

lambda = -80.0;
thresh = 1E-5;

D = p_Data;
%%%%%
%cent = cent(2:end,:);
%%%%%
count = 0;
% while true
% % embed on a circle
% % %measure how close every two consecutive frames:
d = dist2(D,D);
d = [diag(d,1);d(1,end)];
d = d/sum(d);
%accumlate distances so that last value is 1:
accumlate = tril(ones(length(d),length(d)),0);
d = accumlate*d;
d = [0;d(1:end-1)];
%%%
% frames=size(D,1);%angletrain{seq};
% d = ((1:frames)'-1)/frames;
%get embedding on circle.
t=d*2*pi;
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
D = G*CF;
disp 'The reconstruction error is:'
err = norm(D-p_Data)
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
