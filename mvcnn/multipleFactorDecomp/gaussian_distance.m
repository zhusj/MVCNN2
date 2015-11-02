function d = gaussian_distance(x,t)

%d = size(t,1)
%n = size(x,1)
% if size(x,2)>1
%     r = sqrt(pdist2(x,t,'cosine'));
% else
%     r = abs(bsxfun(@minus,x,t'));
% end
thresh = 1E-15;
r = [];
for i = 1:size(t,1)
	r = [r;bsxfun(@minus,x,t(i,:))];
end
r = r';
Q = r*r';
[U,S,V] = svd(Q);
s = diag(S);
l = min(find(s>thresh))
if isempty(l)
    W = V(:,end-d+1:end);
else
    W = V(:,l-d:l-1);
end
% if d==n
%     W = eye(d);
% else
%     
%     
%     Q = r*r';
%     [U,S,V] = svd(Q);
%     s = diag(S);
%     l = min(find(s<thresh))
%     if isempty(l)
%         W = V(:,end-d+1:end);
%     else
%         W = V(:,l-d:l);
%     end
% end

%M = W'*W;
%A = r'*M*r;
m = r'*W;
m = m.^2;
d = sum(m,2);
d = reshape(d,size(x,1),size(t,1));
d = exp(d);



