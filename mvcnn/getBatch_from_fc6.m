function [im, batch_labels] = getBatch_from_fc6(batch,data,labels)

im(1,1,:,:) = single(data(batch,:)');
batch_labels = labels(batch);