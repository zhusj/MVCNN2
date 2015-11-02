
function labels = load_class_labels(GCF,class_dim_nm)
%LEARN_CLASS_LABEL Summary of this function goes here
%   Detailed explanation goes here
labels = {GCF.(class_dim_nm)};
labels = cell2mat(labels);
labels = labels(:);
end
