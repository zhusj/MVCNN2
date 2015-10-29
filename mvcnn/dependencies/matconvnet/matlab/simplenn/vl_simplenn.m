function res = vl_simplenn(net, x, dzdy, res, varargin)
% VL_SIMPLENN  Evaluates a simple CNN
%   RES = VL_SIMPLENN(NET, X) evaluates the convnet NET on data X.
%   RES = VL_SIMPLENN(NET, X, DZDY) evaluates the convnent NET and its
%   derivative on data X and output derivative DZDY.
%
%   The network has a simple (linear) topology, i.e. the computational
%   blocks are arranged in a sequence of layers. Please note that
%   there is no need to use this wrapper, which is provided for
%   convenience. Instead, the individual CNN computational blocks can
%   be evaluated directly, making it possible to create significantly
%   more complex topologies, and in general allowing greater
%   flexibility.
%
%   The NET structure contains two fields:
%
%   - net.layers: the CNN layers.
%   - net.normalization: information on how to normalize input data.
%
%   The network expects the data X to be already normalized. This
%   usually involves rescaling the input image(s) and subtracting a
%   mean.
%
%   RES is a structure array with one element per network layer plus
%   one representing the input. So RES(1) refers to the zeroth-layer
%   (input), RES(2) refers to the first layer, etc. Each entry has
%   fields:
%
%   - res(i+1).x: the output of layer i. Hence res(1).x is the network
%     input.
%
%   - res(i+1).aux: auxiliary output data of layer i. For example,
%     dropout uses this field to store the dropout mask.
%
%   - res(i+1).dzdx: the derivative of the network output relative to
%     variable res(i+1).x, i.e. the output of layer i. In particular
%     res(1).dzdx is the derivative of the network output with respect
%     to the network input.
%
%   - res(i+1).dzdw: the derivative of the network output relative to
%     the parameters of layer i. It can be a cell array for multiple
%     parameters.
%
%   net.layers is a cell array of network layers. The following
%   layers, encapsulating corresponding functions in the toolbox, are
%   supported:
%
%   Convolutional layer::
%     The convolutional layer wraps VL_NNCONV(). It has fields:
%
%     - layer.type = 'conv'
%     - layer.weights = {filters, biases}
%     - layer.stride: the sampling stride (usually 1).
%     - layer.pad: the padding (usually 0).
%
%   Convolution transpose layer::
%     The convolution transpose layer wraps VL_NNCONVT(). It has fields:
%
%     - layer.type = 'convt'
%     - layer.weights = {filters, biases}
%     - layer.upsample: the upsampling factor.
%     - layer.crop: the amount of output cropping.
%
%   Max pooling layer::
%     The max pooling layer wraps VL_NNPOOL(). It has fields:
%
%     - layer.type = 'pool'
%     - layer.method: pooling method ('max' or 'avg').
%     - layer.pool: the pooling size.
%     - layer.stride: the sampling stride (usually 1).
%     - layer.pad: the padding (usually 0).
%
%   Normalization layer::
%     The normalization layer wraps VL_NNNORMALIZE(). It has fields
%
%     - layer.type = 'normalize'
%     - layer.param: the normalization parameters.
%
%   Spatial normalization layer::
%     This is similar to the layer above, but wraps VL_NNSPNORM():
%
%     - layer.type = 'spnorm'
%     - layer.param: the normalization parameters.
%
%   Batch normalization layer::
%     This layer wraps VL_NNBNORM(). It has fields:
%
%     - layer.type = 'bnorm'
%     - layer.weights = {multipliers, biases}.
%
%   ReLU and Sigmoid layers::
%     The ReLU layer wraps VL_NNRELU(). It has fields:
%
%     - layer.type = 'relu'
%
%     The sigmoid layer is the same, but for the sigmoid function, with
%     `relu` replaced by `sigmoid`.
%
%   Dropout layer::
%     The dropout layer wraps VL_NNDROPOUT(). It has fields:
%
%     - layer.type = 'dropout'
%     - layer.rate: the dropout rate.
%
%   Softmax layer::
%     The softmax layer wraps VL_NNSOFTMAX(). It has fields
%
%     - layer.type = 'softmax'
%
%   Log-loss layer::
%     The log-loss layer wraps VL_NNLOSS(). It has fields:
%
%     - layer.type = 'loss'
%     - layer.class: the ground-truth class.
%
%   Softmax-log-loss layer::
%     The softmax-log-loss layer wraps VL_NNSOFTMAXLOSS(). It has
%     fields:
%
%     - layer.type = 'softmaxloss'
%     - layer.class: the ground-truth class.
%
%   P-dist layer::
%     The pdist layer wraps VL_NNPDIST(). It has fields:
%
%     - layer.type = 'pdist'
%     - layer.p = P parameter of the P-distance
%     - layer.noRoot = whether to raise the distance to the P-th power
%     - layer.epsilon = regularization parameter for the derivatives
%
%   Custom layer::
%     This can be used to specify custom layers.
%
%     - layer.type = 'custom'
%     - layer.forward: a function handle computing the block.
%     - layer.backward: a function handle computing the block derivative.
%
%     The first function is called as res(i+1) = forward(layer, res(i), res(i+1))
%     where res() is the struct array specified before. The second function is
%     called as res(i) = backward(layer, res(i), res(i+1)). Note that the
%     `layer` structure can contain additional fields if needed.

% Copyright (C) 2014 Andrea Vedaldi.
% All rights reserved.
%
% This file is part of the VLFeat library and is made available under
% the terms of the BSD license (see the COPYING file).

opts.res = [] ;
opts.conserveMemory = false ;
opts.sync = false ;
opts.disableDropout = false ;
opts.freezeDropout = false ;
opts.accumulate = false ;
opts.cudnn = true ;
opts.backPropDepth = +inf ;
opts.addSupervision = false;
opts.views = 12;

opts = vl_argparse(opts, varargin);

n = numel(net.layers) ;

if (nargin <= 2) || isempty(dzdy)
  doder = false ;
else
  doder = true ;
end

if opts.cudnn
  cudnn = {'CuDNN'} ;
else
  cudnn = {'NoCuDNN'} ;
end

gpuMode = isa(x, 'gpuArray') ;

if nargin <= 3 || isempty(res)
  res = struct(...
    'x', cell(1,n+1), ...
    'dzdx', cell(1,n+1), ...
    'dzdw', cell(1,n+1), ...
    'aux', cell(1,n+1), ...
    'time', num2cell(zeros(1,n+1)), ...
    'backwardTime', num2cell(zeros(1,n+1)), ...
    'pose', cell(1,n+1)) ;
end
res(1).x = x ;


%%%%%%%%%%%%%%%%
% if opts.addSupervision,
%     res_n = struct(...
%         'x', cell(1,2), ...
%         'dzdx', cell(1,2), ...
%         'dzdw', cell(1,2), ...
%         'aux', cell(1,2), ...
%         'time', num2cell(zeros(1,2)), ...
%         'backwardTime', num2cell(zeros(1,2)), ...
%         'pose', cell(1,2)) ;
%     opts.scale = 1;
%     opts.weightDecay = 1;
% %     initBias= 0.1;
%     % net = add_block(net, opts, 9, 1, 1, 4096, 60, 1, 0, 0.1); 
% %     net_n.layers{1} = struct('type', 'conv', 'name', 'fc_p', ...
% %                                'weights', {{0.01/opts.scale * randn(1, 1, 256, 60, 'single'), ...
% %                                initBias*ones(1,60,'single')}}, ...
% %                                'stride', 1, ...
% %                                'pad', 0, ...
% %                                'learningRate', [1 2], ...
% %                                'weightDecay', [opts.weightDecay 0]) ;
%     % net.layers{end+1} = struct('type', 'relu', 'name', 'relu_p') ;
%     net_n.layers{1} = struct('type', 'softmaxloss', 'name', 'loss') ;
%     pose = 1:12;
%     poses = repmat(pose,1,size(x,4)/12);
%     net_n.layers{1}.class = poses;
%     net_n = vl_simplenn_move(net_n, 'gpu') ;
% end
if opts.addSupervision && size(x,4)~=1
%       res = struct(...
%     'x', cell(1,n+1), ...
%     'dzdx', cell(1,n+1), ...
%     'dzdw', cell(1,n+1), ...
%     'aux', cell(1,n+1), ...
%     'time', num2cell(zeros(1,n+1)), ...
%     'backwardTime', num2cell(zeros(1,n+1)), ...
%     'pose', cell(1,n+1)) ;
    res_n = struct(...
        'x', cell(1,3), ...
        'dzdx', cell(1,3), ...
        'dzdw', cell(1,3), ...
        'aux', cell(1,3), ...
        'time', num2cell(zeros(1,3)), ...
        'backwardTime', num2cell(zeros(1,3)), ...
        'pose', cell(1,3)) ;
%     res_n = struct(...
%         'x', cell(size(x,4)/12,3), ...
%         'dzdx', cell(size(x,4)/12,3), ...
%         'dzdw', cell(size(x,4)/12,3), ...
%         'aux', cell(size(x,4)/12,3), ...
%         'time', num2cell(zeros(size(x,4)/12,3)), ...
%         'backwardTime', num2cell(zeros(size(x,4)/12,3)), ...
%         'pose', cell(size(x,4)/12,3)) ;
    opts.scale = 1;
    opts.weightDecay = 1;
    initBias= 0.1;
    % net = add_block(net, opts, 9, 1, 1, 4096, 60, 1, 0, 0.1); 
    net_n.layers{1} = struct('type', 'conv', 'name', 'fc_p', ...
                               'weights', {{0.01/opts.scale * randn(1, 1, 512, opts.views, 'single'), ...
                               initBias*ones(1,opts.views,'single')}}, ...
                               'stride', 1, ...
                               'pad', 0, ...
                               'learningRate', [10 20], ...
                               'weightDecay', [opts.weightDecay 0]) ;
    % net.layers{end+1} = struct('type', 'relu', 'name', 'relu_p') ;
    net_n.layers{2} = struct('type', 'softmaxloss', 'name', 'loss') ;
    pose = 1:opts.views;
    poses = repmat(pose,1,size(x,4)/opts.views);
    net_n.layers{2}.class = poses;
    net_n = vl_simplenn_move(net_n, 'gpu') ;
end
%%%%%%%%%%%%%%%%



for i=1:n
  l = net.layers{i} ;
  res(i).time = tic ;
  switch l.type
    case 'conv'
      if isfield(l, 'weights')
%         if i == 25
%             res(i+1).x = vl_nnconv(res(14).x, l.weights{1}, l.weights{2}, ...
%                        'pad', l.pad, 'stride', l.stride, ...
%                        cudnn{:}) ;
%         else
%             res(i+1).x = vl_nnconv(res(i).x, l.weights{1}, l.weights{2}, ...
%                        'pad', l.pad, 'stride', l.stride, ...
%                        cudnn{:}) ;
%         end
        res(i+1).x = vl_nnconv(res(i).x, l.weights{1}, l.weights{2}, ...
                               'pad', l.pad, 'stride', l.stride, ...
                               cudnn{:}) ;
      else
           if i == 23
          %%%%% L2 norm %%%%%%%%%%%%%%%%%%
              res_23_sum = sqrt(sum(sum(sum(sum(abs(res(23).x).^2)))));
              res_20_sum = sqrt(sum(sum(sum(sum(abs(res(20).x).^2)))));
              res(i+1).x = vl_nnconv(res(20).x/res_20_sum*10 + res(23).x/res_23_sum*10, l.filters, l.biases, ...
                 'pad', l.pad, 'stride', l.stride, ...
                 cudnn{:}) ;
           else
              res(i+1).x = vl_nnconv(res(i).x, l.filters, l.biases, ...
                       'pad', l.pad, 'stride', l.stride, ...
                       cudnn{:}) ;
           end
%           res(i+1).x = vl_nnconv(res(i).x, l.filters, l.biases, ...
%                        'pad', l.pad, 'stride', l.stride, ...
%                        cudnn{:}) ;

%         if isfield(l,'branch')
%             if strcmp(l.branch, 'pose')
%                 res(i+1).x = vl_nnconv(res(i-2).x, l.filters, l.biases, ...
%                            'pad', l.pad, 'stride', l.stride, ...
%                            cudnn{:}) ;
%             
%             else
%                 res(i+1).x = vl_nnconv(res(i).x, l.filters, l.biases, ...
%                                        'pad', l.pad, 'stride', l.stride, ...
%                                        cudnn{:}) ;
%             end
%         end
      end
    case 'convt'
      if isfield(l, 'weights')
        res(i+1).x = vl_nnconvt(res(i).x, l.weights{1}, l.weights{2}, ...
                               'crop', l.crop, 'upsample', l.upsample, ...
                               'numGroups', l.numGroups, cudnn{:}) ;
      else
        res(i+1).x = vl_nnconvt(res(i).x, l.filters, l.biases, ...
                               'crop', l.pad, 'upsample', l.upsample, ...
                               'numGroups', l.numGroups, cudnn{:}) ;
      end
    case 'pool'
      res(i+1).x = vl_nnpool(res(i).x, l.pool, ...
                             'pad', l.pad, 'stride', l.stride, ...
                             'method', l.method, ...
                             cudnn{:}) ;
    case 'normalize'
      res(i+1).x = vl_nnnormalize(res(i).x, l.param) ;
    case 'softmax'
      res(i+1).x = vl_nnsoftmax(res(i).x) ;
    case 'loss'
      res(i+1).x = vl_nnloss(res(i).x, l.class) ;
    case 'softmaxloss'
%       lambda_1 = 0.5;
%       lambda_2 = 0.5;
%       angles = [30,60,90,120,150,180,210,240,270,300,330,360]';
%       error = zeros(1,5);
%       x = gather(res(15).pose);
%       for j = 1:5
%           error(j) = manifold_alignment(x(:,:,j)', angles, 'correlation'); 
%       end
%       res(i+1).x = lambda_1 * vl_nnsoftmaxloss(res(i).x, l.class) + lambda_2 * mean(error) ;
      res(i+1).x = vl_nnsoftmaxloss(res(i).x, l.class) ;
    case 'relu'
      if isfield(l, 'leak'), leak = {'leak', l.leak} ; else leak = {} ; end
%       leak = {'leak', 0.01};
      res(i+1).x = vl_nnrelu(res(i).x,[],leak{:}) ;
    case 'sigmoid'
      res(i+1).x = vl_nnsigmoid(res(i).x) ;
    case 'noffset'
      res(i+1).x = vl_nnnoffset(res(i).x, l.param) ;
    case 'spnorm'
      res(i+1).x = vl_nnspnorm(res(i).x, l.param) ;
    case 'dropout'
      if opts.disableDropout
        res(i+1).x = res(i).x ;
      elseif opts.freezeDropout
        [res(i+1).x, res(i+1).aux] = vl_nndropout(res(i).x, 'rate', l.rate, 'mask', res(i+1).aux) ;
      else
        [res(i+1).x, res(i+1).aux] = vl_nndropout(res(i).x, 'rate', l.rate) ;
      end
    case 'bnorm'
      if isfield(l, 'weights')
        res(i+1).x = vl_nnbnorm(res(i).x, l.weights{1}, l.weights{2}) ;
      else
        res(i+1).x = vl_nnbnorm(res(i).x, l.filters, l.biases) ;
      end
    case 'pdist'
      res(i+1) = vl_nnpdist(res(i).x, l.p, 'noRoot', l.noRoot, 'epsilon', l.epsilon) ;
    case 'custom'
%         res(i+1) = l.forward(l, res(14), res(i+1)) ;
      res(i+1) = l.forward(l, res(i), res(i+1)) ;
    otherwise
      error('Unknown layer type %s', l.type) ;
  end
  % optionally forget intermediate results
  forget = opts.conserveMemory ;
  forget = forget & (~doder || strcmp(l.type, 'relu')) ;
  forget = forget & ~(strcmp(l.type, 'loss') || strcmp(l.type, 'softmaxloss')) ;
  forget = forget & (~isfield(l, 'rememberOutput') || ~l.rememberOutput) ;
  if forget && i~=14
    res(i).x = [] ;
  end
  if gpuMode & opts.sync
    % This should make things slower, but on MATLAB 2014a it is necessary
    % for any decent performance.
    wait(gpuDevice) ;
  end
  res(i).time = toc(res(i).time) ;
end

%%%%%%%%%%%%%%%%%%%%%%%%
if doder && opts.addSupervision
    for i=1:2
      l = net_n.layers{i} ;
%       for j=1:5
          res_n(i).time = tic ;
          res_n(1).x = res(14).x;
%           res_n(j,1).x = res(14).x(:,:,:,(j-1)*12+1:j*12);
          switch l.type
            case 'conv'
              if isfield(l, 'weights')
                res_n(i+1).x = vl_nnconv(res_n(i).x, l.weights{1}, l.weights{2}, ...
                                       'pad', l.pad, 'stride', l.stride, ...
                                       cudnn{:}) ;
              else
                  res(i+1).x = vl_nnconv(res_n(i).x, l.filters, l.biases, ...
                               'pad', l.pad, 'stride', l.stride, ...
                               cudnn{:}) ;
              end
            case 'softmaxloss'
              res_n(i+1).x = vl_nnsoftmaxloss(res_n(i).x, l.class) ;
            otherwise
              error('Unknown layer type %s', l.type) ;
          end
          % optionally forget intermediate results
          forget = opts.conserveMemory ;
          forget = forget & (~doder || strcmp(l.type, 'relu')) ;
          forget = forget & ~(strcmp(l.type, 'loss') || strcmp(l.type, 'softmaxloss')) ;
          forget = forget & (~isfield(l, 'rememberOutput') || ~l.rememberOutput) ;
          if forget
            res_n(i).x = [] ;
          end
          if gpuMode & opts.sync
            % This should make things slower, but on MATLAB 2014a it is necessary
            % for any decent performance.
            wait(gpuDevice) ;
          end
          res_n(i).time = toc(res_n(i).time) ;
%       end
    end
end

if doder && opts.addSupervision  
  for i=2:-1:1
    l = net_n.layers{i} ;
%     for j=1:5
        res_n(3).dzdx = dzdy ;
        res_n(i).backwardTime = tic ;
        switch l.type
          case 'conv'
            if ~opts.accumulate
              if isfield(l, 'weights')
                [res_n(i).dzdx, res_n(i).dzdw{1}, res_n(i).dzdw{2}] = ...
                    vl_nnconv(res_n(i).x, l.weights{1}, l.weights{2}, ...
                              res_n(i+1).dzdx, ...
                              'pad', l.pad, 'stride', l.stride, ...
                              cudnn{:}) ;
              else
                % Legacy code: will go
                [res(i).dzdx, res(i).dzdw{1}, res(i).dzdw{2}] = ...
                    vl_nnconv(res(i).x, l.filters, l.biases, ...
                              res(i+1).dzdx, ...
                              'pad', l.pad, 'stride', l.stride, ...
                              cudnn{:}) ;
              end
            else
              dzdw = cell(1,2) ;
              if isfield(l, 'weights')
                [res(i).dzdx, dzdw{1}, dzdw{2}] = ...
                    vl_nnconv(res(i).x, l.weights{1}, l.weights{2}, ...
                              res(i+1).dzdx, ...
                              'pad', l.pad, 'stride', l.stride, ...
                              cudnn{:}) ;
              else
                % Legacy code: will go
                [res(i).dzdx, dzdw{1}, dzdw{2}] = ...
                    vl_nnconv(res(i).x, l.filters, l.biases, ...
                              res(i+1).dzdx, ...
                              'pad', l.pad, 'stride', l.stride, ...
                              cudnn{:}) ;
              end
              for jj=1:2
                res(i).dzdw{jj} = res(i).dzdw{jj} + dzdw{jj} ;
              end
              clear dzdw ;
            end
          case 'softmaxloss'
            res_n(i).dzdx = vl_nnsoftmaxloss(res_n(i).x, l.class, res_n(i+1).dzdx) ;
        end
        if opts.conserveMemory
          res_n(i+1).dzdx = [] ;
        end
        if gpuMode & opts.sync
          wait(gpuDevice) ;
        end
        res_n(i).backwardTime = toc(res_n(i).backwardTime) ;
%     end
  end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if doder
  res(n+1).dzdx = dzdy ;
  for i=n:-1:max(1, n-opts.backPropDepth+1)
    l = net.layers{i} ;
    res(i).backwardTime = tic ;
    switch l.type
      case 'conv'
        if ~opts.accumulate
          if isfield(l, 'weights')
%               if i == 13
%                   [res(i).dzdx, res(i).dzdw{1}, res(i).dzdw{2}] = ...
%                     vl_nnconv(res(i).x, l.weights{1}, l.weights{2}, ...
%                               res(i+1).dzdx + res(25).dzdx, ...
%                               'pad', l.pad, 'stride', l.stride, ...
%                               cudnn{:}) ;              
%               elseif i == 25
%                   [res_p.dzdx, res(i).dzdw{1}, res(i).dzdw{2}] = ...
%                     vl_nnconv(res(14).x, l.weights{1}, l.weights{2}, ...
%                               res(i+1).dzdx, ...
%                               'pad', l.pad, 'stride', l.stride, ...
%                               cudnn{:}) ;
%               else
%                   [res(i).dzdx, res(i).dzdw{1}, res(i).dzdw{2}] = ...
%                     vl_nnconv(res(i).x, l.weights{1}, l.weights{2}, ...
%                               res(i+1).dzdx, ...
%                               'pad', l.pad, 'stride', l.stride, ...
%                               cudnn{:}) ;
%               end
%             if i == 14
%                 [res_p.dzdx, res(i).dzdw{1}, res(i).dzdw{2}] = ...
%                   vl_nnconv(res(i).x, l.weights{1}, l.weights{2}, ...
%                             res(i+1).dzdx, ...
%                             'pad', l.pad, 'stride', l.stride, ...
%                             cudnn{:}) ;
%             else
%                 [res(i).dzdx, res(i).dzdw{1}, res(i).dzdw{2}] = ...
%                   vl_nnconv(res(i).x, l.weights{1}, l.weights{2}, ...
%                             res(i+1).dzdx, ...
%                             'pad', l.pad, 'stride', l.stride, ...
%                             cudnn{:}) ;
%             end

            [res(i).dzdx, res(i).dzdw{1}, res(i).dzdw{2}] = ...
                vl_nnconv(res(i).x, l.weights{1}, l.weights{2}, ...
                          res(i+1).dzdx, ...
                          'pad', l.pad, 'stride', l.stride, ...
                          cudnn{:}) ;
          else
            % Legacy code: will go
            if i == 13
                res_14_dzdx_sum = sum(sum(sum(sum(abs(res(14).dzdx)))));
                pose_dzdx_sum = sum(sum(sum(sum(abs(res_n(1).dzdx)))));

            %%%%% L2 norm %%%%%%%%%%%%%%%%%%
%                 res_14_dzdx_sum = sqrt(sum(sum(sum(sum(abs(res(14).dzdx).^2)))));
%                 pose_dzdx_sum = sqrt(sum(sum(sum(sum(abs(res_n(1).dzdx).^2)))));

                [res(i).dzdx, res(i).dzdw{1}, res(i).dzdw{2}] = ...
                    vl_nnconv(res(i).x, l.filters, l.biases, ...
                              res(14).dzdx + res_n(1).dzdx/pose_dzdx_sum*res_14_dzdx_sum, ...
                              'pad', l.pad, 'stride', l.stride, ...
                              cudnn{:}) ;
       %%%%%%%%%% Decay
%                 [res(i).dzdx, res(i).dzdw{1}, res(i).dzdw{2}] = ...
%                     vl_nnconv(res(i).x, l.filters, l.biases, ...
%                               res(14).dzdx + 1*(1-net.epoch/net.numEpochs)*res_n(1).dzdx/pose_dzdx_sum*res_14_dzdx_sum, ...
%                               'pad', l.pad, 'stride', l.stride, ...
%                               cudnn{:}) ;
            else
                [res(i).dzdx, res(i).dzdw{1}, res(i).dzdw{2}] = ...
                    vl_nnconv(res(i).x, l.filters, l.biases, ...
                              res(i+1).dzdx, ...
                              'pad', l.pad, 'stride', l.stride, ...
                              cudnn{:}) ;
            end

%             [res(i).dzdx, res(i).dzdw{1}, res(i).dzdw{2}] = ...
%                 vl_nnconv(res(i).x, l.filters, l.biases, ...
%                           res(i+1).dzdx, ...
%                           'pad', l.pad, 'stride', l.stride, ...
%                           cudnn{:}) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%% pose manifold error %%%%%%%%%%%%%%%%%%%%%%%5
%             if i == 13
% %               angles = linspace(0,2*pi,12)';
%               [sz1, sz2, sz3, sz4] = size(res(i+1).dzdx);
%               theta = linspace(0,2*pi,12)';
%               pts_unit_circle = [cos(theta), sin(theta)];
%               pose = gather(res(15).pose);
%               pose_dzdx = [];
%               
%               
%               for j = 1:sz4/12
%                   x = pose(:,:,1)';
%                   %using pca (from drtoolbox)
%                   [x_proj,mapping]= pca2(x,2);
%                   x_proj_n = (x_proj - repmat(mean(x_proj),...
%                         [size(x_proj,1) 1])) ./ repmat(std(x_proj),[size(x_proj,1) 1]);
%         
% %                   x_proj_n = normalizePoints(x_proj);
%                   pose_gradient = 2*(x_proj_n-pts_unit_circle);
%                   
%                   pose_gradient_n = (pose_gradient + repmat(mean(x_proj),...
%                         [size(pose_gradient,1) 1])) ./ repmat(std(x_proj),[size(pose_gradient,1) 1]);
%                   % back projection
%                   try
%                       pose_gradient_backproj = (pinv(mapping.M)' * pose_gradient_n')';
%                   catch
%                       flag=1;
%                   end
% %                   pose_gradient_backproj = bsxfun(@plus, pose_gradient_backproj, mapping.mean);
%                   
% %                   pose_gradient_backproj_n = pose_gradient_backproj/sum(sum(sum(sum(abs(pose_gradient_backproj)))));
% 
%                 
%                   pose_gradient_backproj_r = reshape(pose_gradient_backproj',sz1,sz2,sz3,12);
%                   
%                   pose_dzdx = cat(4,pose_dzdx, pose_gradient_backproj_r);
%               end
%               
% %               pose_error = zeros(1,5);
% %               x = gather(res(15).pose);
% %               for j = 1:5
% %                   pose_error(j) = manifold_alignment(x(:,:,j)', angles, 'correlation'); 
% %               end
% %               [sz1, sz2, sz3, sz4] = size(res(i+1).dzdx);
% %               pose_error = reshape(repmat(pose_error,12,1),60,1);
% %               pose_gradient = permute(repmat(pose_error,1,sz1,sz2,sz3),[2,3,4,1])/(sz1*sz2*sz3);
% %               [res(i).dzdx, res(i).dzdw{1}, res(i).dzdw{2}] = ...
% %                 vl_nnconv(res(i).x, l.filters, l.biases, ...
% %                           res(i+1).dzdx, ...
% %                           'pad', l.pad, 'stride', l.stride, ...
% %                           cudnn{:}) ;
%               res_14_dzdx_sum = sum(sum(sum(sum(abs(res(i+1).dzdx)))));
%               pose_dzdx_sum = sum(sum(sum(sum(abs(pose_dzdx)))));
%               pose_dzdx = pose_dzdx/pose_dzdx_sum * res_14_dzdx_sum/pose_dzdx_sum;
%               try
%                   [res(i).dzdx, res(i).dzdw{1}, res(i).dzdw{2}] = ...
%                     vl_nnconv(res(i).x, l.filters, l.biases, ...
%                               res(i+1).dzdx + pose_dzdx, ...
%                               'pad', l.pad, 'stride', l.stride, ...
%                               cudnn{:}) ;
%               catch
%                   flag =1;
%               end
%             else
%               [res(i).dzdx, res(i).dzdw{1}, res(i).dzdw{2}] = ...
%                 vl_nnconv(res(i).x, l.filters, l.biases, ...
%                           res(i+1).dzdx, ...
%                           'pad', l.pad, 'stride', l.stride, ...
%                           cudnn{:}) ;
%             end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          end
        else
          dzdw = cell(1,2) ;
          if isfield(l, 'weights')
            [res(i).dzdx, dzdw{1}, dzdw{2}] = ...
                vl_nnconv(res(i).x, l.weights{1}, l.weights{2}, ...
                          res(i+1).dzdx, ...
                          'pad', l.pad, 'stride', l.stride, ...
                          cudnn{:}) ;
          else
            % Legacy code: will go
            [res(i).dzdx, dzdw{1}, dzdw{2}] = ...
                vl_nnconv(res(i).x, l.filters, l.biases, ...
                          res(i+1).dzdx, ...
                          'pad', l.pad, 'stride', l.stride, ...
                          cudnn{:}) ;
          end
%             if isfield(l,'branch')
%                 if strcmp(l.branch, 'pose')
%                     [res(i).dzdx, dzdw{1}, dzdw{2}] = ...
%                         vl_nnconv(res(i).x, l.filters, l.biases, ...
%                                   res(i+1).dzdx, ...
%                                   'pad', l.pad, 'stride', l.stride, ...
%                                   cudnn{:}) ;              
%                 else
%                 [res(i).dzdx, dzdw{1}, dzdw{2}] = ...
%                     vl_nnconv(res(i).x, l.filters, l.biases, ...
%                               res(i+1).dzdx, ...
%                               'pad', l.pad, 'stride', l.stride, ...
%                               cudnn{:}) ;
%                 end
%             end

          for j=1:2
            res(i).dzdw{j} = res(i).dzdw{j} + dzdw{j} ;
          end
          clear dzdw ;
        end

      case 'convt'
        if ~opts.accumulate
          if isfield(l, 'weights')
            [res(i).dzdx, res(i).dzdw{1}, res(i).dzdw{2}] = ...
                vl_nnconvt(res(i).x, l.weights{1}, l.weights{2}, ...
                          res(i+1).dzdx, ...
                          'crop', l.crop, 'upsample', l.upsample, ...
                          'numGroups', l.numGroups, cudnn{:}) ;
          else
            % Legacy code: will go
            [res(i).dzdx, res(i).dzdw{1}, res(i).dzdw{2}] = ...
                vl_nnconvt(res(i).x, l.filters, l.biases, ...
                         res(i+1).dzdx, ...
                          'crop', l.crop, 'upsample', l.upsample, ...
                          'numGroups', l.numGroups, cudnn{:}) ;          end
        else
          dzdw = cell(1,2) ;
          if isfield(l, 'weights')
            [res(i).dzdx, dzdw{1}, dzdw{2}] = ...
                vl_nnconvt(res(i).x, l.weights{1}, l.weights{2}, ...
                          res(i+1).dzdx, ...
                          'crop', l.crop, 'upsample', l.upsample, ...
                          'numGroups', l.numGroups, cudnn{:}) ;          else
            % Legacy code: will go
            [res(i).dzdx, dzdw{1}, dzdw{2}] = ...
                vl_nnconvt(res(i).x, l.filters, l.biases, ...
                          res(i+1).dzdx, ...
                          'crop', l.crop, 'upsample', l.upsample, ...
                          'numGroups', l.numGroups, cudnn{:}) ;          end
          for j=1:2
            res(i).dzdw{j} = res(i).dzdw{j} + dzdw{j} ;
          end
          clear dzdw ;
        end

      case 'pool'
        res(i).dzdx = vl_nnpool(res(i).x, l.pool, res(i+1).dzdx, ...
                                'pad', l.pad, 'stride', l.stride, ...
                                'method', l.method, ...
                                cudnn{:}) ;
      case 'normalize'
        res(i).dzdx = vl_nnnormalize(res(i).x, l.param, res(i+1).dzdx) ;
      case 'softmax'
        res(i).dzdx = vl_nnsoftmax(res(i).x, res(i+1).dzdx) ;
      case 'loss'
        res(i).dzdx = vl_nnloss(res(i).x, l.class, res(i+1).dzdx) ;
      case 'softmaxloss'
        %%%%%% my change to accomendate two softmaxloss
%         res(i+1).dzdx = dzdy ;
        %%%%%
        res(i).dzdx = vl_nnsoftmaxloss(res(i).x, l.class, res(i+1).dzdx) ;
      case 'relu'
        if isfield(l, 'leak'), leak = {'leak', l.leak} ; else leak = {} ; end
%         leak = {'leak', 0.01};
        if ~isempty(res(i).x)
          res(i).dzdx = vl_nnrelu(res(i).x, res(i+1).dzdx, leak{:}) ;
        else
          % if res(i).x is empty, it has been optimized away, so we use this
          % hack (which works only for ReLU):
          res(i).dzdx = vl_nnrelu(res(i+1).x, res(i+1).dzdx, leak{:}) ;
        end
      case 'sigmoid'
        res(i).dzdx = vl_nnsigmoid(res(i).x, res(i+1).dzdx) ;
      case 'noffset'
        res(i).dzdx = vl_nnnoffset(res(i).x, l.param, res(i+1).dzdx) ;
      case 'spnorm'
        res(i).dzdx = vl_nnspnorm(res(i).x, l.param, res(i+1).dzdx) ;
      case 'dropout'
        if opts.disableDropout
          res(i).dzdx = res(i+1).dzdx ;
        else
          res(i).dzdx = vl_nndropout(res(i).x, res(i+1).dzdx, ...
                                     'mask', res(i+1).aux) ;
        end
      case 'bnorm'
        if ~opts.accumulate
          if isfield(l, 'weights')
            [res(i).dzdx, res(i).dzdw{1}, res(i).dzdw{2}] = ...
                vl_nnbnorm(res(i).x, l.weights{1}, l.weights{2}, ...
                           res(i+1).dzdx) ;
          else
            [res(i).dzdx, res(i).dzdw{1}, res(i).dzdw{2}] = ...
                vl_nnbnorm(res(i).x, l.filters, l.biases, ...
                           res(i+1).dzdx) ;
          end
        else
          dzdw = cell(1,2) ;
          if isfield(l, 'weights')
            [res(i).dzdx, dzdw{1}, dzdw{2}] = ...
                vl_nnbnorm(res(i).x, l.weights{1}, l.weights{2}, ...
                           res(i+1).dzdx) ;
          else
            [res(i).dzdx, dzdw{1}, dzdw{2}] = ...
                vl_nnbnorm(res(i).x, l.filters, l.biases, ...
                           res(i+1).dzdx) ;
          end
          for j=1:2
            res(i).dzdw{j} = res(i).dzdw{j} + dzdw{j} ;
          end
          clear dzdw ;
        end
      case 'pdist'
        res(i).dzdx = vl_nnpdist(res(i).x, l.p, res(i+1).dzdx, ...
                                 'noRoot', l.noRoot, 'epsilon', l.epsilon) ;
      case 'custom'
        res(i) = l.backward(l, res(i), res(i+1)) ;
    end
    if opts.conserveMemory
      res(i+1).dzdx = [] ;
    end
    if gpuMode & opts.sync
      wait(gpuDevice) ;
    end
    res(i).backwardTime = toc(res(i).backwardTime) ;
  end
end






