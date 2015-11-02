function L = get_matrix_label(lbl,c)
    if ~exist('c','var') || c==0
        c = max(lbl);%length(unique(lbl));
    end
    len = length(lbl);
    L = zeros(len,c);
    L(sub2ind(size(L),(1:len)',lbl)) = 1;
%     s = sum(L);
%     L(:,~s)=[];
    %L = L(:,1:end-1);
    
    