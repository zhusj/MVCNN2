function [cls,err_final]=classify_vector(samples,bases, labels)

S = cosine_similarity(samples, bases);
[~,m] = max(S,[],2);
cls = labels(m);
return;
% if(p_use_wts)
%     [ltr,D] = knnsearch(letter_bases,l','K',5,'Distance','seuclidean','Scale',ltr_wts);
%     %ltr = ceil(ltr/mltr)%????????????
%     ltr = mod(ltr-1,numofletters)+1;
%     D
%     ltr_unq = unique(ltr);%,'stable');
%     ltr_rep = histc(ltr,ltr_unq)
%     [dmy,idx] = max(ltr_rep)%ToDo if more than one max value then should choose based on sum of wights, currently chooose the lowest ranked one
%     letter_est_results(li,si) = ltr_unq(idx);
% else
    %cls = knnclassify(v,bases,labels,1,'cosine');
    %cls = mod(cls-1,26)+1;
    
    %ltr = mod(best_letter-1,train_seq_dim(1))+1;%knnclassify(l',letter_bases,ltr_labels,1);
    labels = reshape(labels,1,[]);
    cls = zeros(size(samples,1),1);
    err_final = -1*ones(size(samples,1),1);
    for k = 1:size(samples,1)
         cls = knnclassify(samples(:,1:10),bases(:,1:10),labels);%,1,'cosine');
         %cls = mod(cls-1,26)+1;
        %err_final = 0;
         continue;
        err_min = inf;
        v = samples(k,:);
        v = v';
        for i = unique(labels)
            indx = (labels == i);
            A = bases(indx,:)';
            P = A*inv(A'*A)*A';
            pv = P*v;
            err = norm(v - pv);
            d = min(sqrt(dist2(v',A')));
            err = err + d;
            if err<err_min
                cls(k) = i;
                err_min = err;
            end
        end
        err_final(k) = err_min;
    end
end