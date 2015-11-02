function res = branch12Project(input, net, res)
[sz1,sz2,sz3,sz4] = size(input);
for j = 1:sz4/12
    for i =1:12
        x = input(:,:,:,i+(j-1)*12);
        l = net.layers{i+17} ;
        if isfield(l, 'weights')
            out(:,:,:,i) = vl_nnconv(x, l.weights{1}, l.weights{2}, ...
                           'pad', l.pad, 'stride', l.stride, ...
                           'CuDNN') ;
        else
            out(:,:,:,i) = vl_nnconv(x, l.filters, l.biases, ...
                           'pad', l.pad, 'stride', l.stride, ...
                           'CuDNN') ;
        end
    end
    tmp = mean(out,4);
%     tmp = max(out,[],4);
    res(30).x(:,:,:,j) = single(gather(tmp));    
end
res(30).x = gpuArray(res(30).x);
res(30).out = single(res(30).x);