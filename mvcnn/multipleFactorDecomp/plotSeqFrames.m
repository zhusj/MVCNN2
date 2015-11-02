function plotSeqFrames(seq)
m = 60;n = 80;
seqSz = double(size(seq,2));
gridSz = ceil(sqrt(seqSz));
figure;
colormap(gray);
for i = 1: seqSz
    subplot(gridSz,gridSz,i);
    c = reshape(seq(:,i),m,n);
    %c = c/max(max(c));
    imagesc(c);
end