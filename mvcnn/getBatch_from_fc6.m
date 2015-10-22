function [im, batch_labels] = getBatch_from_fc6(batch,data,labels)

im = data(batch,:);
batch_labels = labels(batch);