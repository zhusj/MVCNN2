% learn GRBF mapping
% input:
% D : original data matrix Nxd  N: is the number of samples, d: is the dimensionality of the input images
% X : embedding space representation. Nxe matrix N: number of samples, e: the dimesnionality of the embedding space
% cent: centers: exm e: the dimesnionality of the embedding space

function map=learnmapping_input_output(p_output,p_input,p_labels)
%I have input points, output points and need regression matrix.

%compute centers as the average of latent points for every class
unique_lbls = unique(p_labels);
cent = zeros(numel(unique_lbls),size(p_input,2));
c = 1;
for i = unique_lbls'
    indx = (p_labels==i);
    l = p_input(indx,:);
    cent(c,:) = mean(l);
    c = c+1;
end

lambda = -80.0;
thresh = 1E-5;

D = p_output;
X = p_input;

%implementation of higher dimensional regularization
%cf = (G'G+lambdag)-1 G'y

N=size(D,2);
G =directs_grbf_reg(X,cent);
g =directs_grbf_reg(cent,cent);
% G =directs_grbf_reg(t,icent);
% g =directs_grbf_reg(icent,icent);

F =G'*D;
B =full(G'*G-lambda*g);
map = B\F;

%count = count+1;

% if exist('CF_old','var')
%     diff = norm(CF-CF_old);
%     if diff<thresh || count>20
%         break;
%     end
% end
%CF_old = CF;
D = G*map;
disp 'The reconstruction error is:'
err = norm(D - p_output)
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
