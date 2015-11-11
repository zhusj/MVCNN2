function res = branch12Project(input, net, res, opts)
[sz1,sz2,sz3,sz4] = size(input);
% angles = 0:30:330;
weight = zeros(12,12);
dist = zeros(12,12);
% for i = 1:12
%     for j = 1:12
%         dist(i,j) = exp(-min(abs(angles(i)-angles(j)),360-abs(angles(i)-angles(j)))/30);
%     end
% end
% dist = dist/sum(dist(1,:));

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
    if opts.weighted
    %%%%%%%%%%%%%%%% weighted pooling
    %%%%% to ba added, larger weight for closer views(on top of difference of features)
        tmp1 = squeeze(out);
        for ii = 1:12
            for jj = 1:12
                dist(ii,jj) = norm(tmp1(:,ii)-tmp1(:,jj))/50;
                weight(ii,jj) = exp(-dist(ii,jj));                
            end
        end
        weight_v = sum(weight);
        weight_vector{j} = weight_v/sum(weight_v);
        dist_v = sum(dist);
        dist_vector{j} = dist_v/sum(dist_v);
        tmp2 = tmp1 * weight_vector{j}';
        tmp(1,1,:,:) = tmp2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55


%     tmp1 = tmp1*dist;
%     out(1,1,:,:) = tmp1;
    else
    %%% ave pooling
        tmp = mean(out,4);
    
    %%% max pooling
%       tmp = max(out,[],4);
    end
    res(30).x(:,:,:,j) = single(gather(tmp));    
end

res(30).dist_v = dist_vector;
res(30).x = gpuArray(res(30).x);
res(30).out = single(res(30).x);