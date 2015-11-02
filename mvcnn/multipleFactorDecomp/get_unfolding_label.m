function l = get_unfolding_label(dim,sz)
    l = 1:sz(dim);
    if dim > 1
        l = repmat(l,prod(sz(1:dim-1)),1);
    end
    l = l(:);
    
    if dim < length(sz)
        l = repmat(l,prod(sz(dim+1:end)),1);
    end