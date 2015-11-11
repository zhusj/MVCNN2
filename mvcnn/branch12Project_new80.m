function res = branch12Project_new(input, net, res, opts)
[sz1,sz2,sz3,sz4] = size(input);
% angles = 0:30:330;
% featureWeight = zeros(12,12);
% featureDist = zeros(12,12);
% viewDist = zeros(12,12);
% for i = 1:12
%     for j = 1:12
%         viewDist(i,j) = min(abs(angles(i)-angles(j)),360-abs(angles(i)-angles(j)));
%     end
% end
% viewDist = viewDist/sum(viewDist(1,:));
% viewWeight = exp(-viewDist/0.05);
% viewWeight = viewWeight/sum(viewWeight(1,:));


for j = 1:max(sz4/80,1)
    for i =1:80
        x = input(:,:,:,i+(j-1)*80);
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
        for ii = 1:sz4
            for jj = 1:sz4
                featureDist(ii,jj) = norm(tmp1(:,ii)-tmp1(:,jj))/50;
                featureWeight(ii,jj) = exp(-featureDist(ii,jj));                
            end
        end
        featureWeight = featureWeight/sum(featureWeight(:,1));
        weight_v = sum(featureWeight);%*viewWeight);
        weight_vector{j} = weight_v/sum(weight_v);
        dist_v = sum(featureDist);%viewDist*
        dist_vector{j} = dist_v/sum(dist_v);
        tmp2 = tmp1;% * weight_vector{j}';
        tmp(1,1,:,:) = tmp2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55


%     tmp1 = tmp1*dist;
%     out(1,1,:,:) = tmp1;
    else
    %%% ave pooling
        tmp = mean(out,4);
        
    %%% root mean square pooling
%         tmp1 = squeeze(out);
%         tmp = zeros(2048,1);
%         for i = 1:12
%             tmp = tmp + tmp1(:,i).^2;
%         end
%         tmp = sqrt(tmp)/12;
        
    % concatenate
%         for i =1:12
%             tmp(1,1,500*(i-1)+1:500*i,1) = out(:,:,:,i);
%         end
    %%% max pooling
%       tmp = max(out,[],4);
    end
    if isfield(l, 'weights')
        dim3 = size(l.weights{1},4);
    else
        dim3 = size(l.filters,4);
    end
     res(98).x(1,1,1:dim3,j) = single(gather(tmp));    
end

if opts.weighted
    res(30).dist_v = dist_vector;
end
res(98).x = gpuArray(res(98).x);
res(98).out = single(res(98).x);