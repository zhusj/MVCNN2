function feats = imdb_compute_cnn_features_save( imdbName, model, varargin )
%IMDB_COMPUTE_CNN_FEATURES Compute and save CNN activations features
%
%   imdb:: 'ModelNet40v1'
%       name of a folder under 'data/'
%   model:: 'imagenet-vgg-m'
%       can be either string (model name) or the actual net model
%       model will be searched/saved under 'data/models'
%   `aug`:: 'none'
%       1st field(f|n) indicates whether include flipped copy or not
%       2nd field(s|r) indicates type of region - Square or Rectangle
%       3rd field(1..4) indicates number of levels
%       note: 'none', 'ns1', 'nr1' are equivalent
%   `gpuMode`:: false
%       set to true to compute on GPU
%   `numWorkers`:: 12
%       number of CPU workers, only in use when gpuMode is false
%   `restart`:: false
%       set to true to re-compute all features
%   `readOp`:: @imread_255
%       the operator that reads data from file
%   `normalization`:: true
%       set to false to turn off all forms of post-processing
%       (normalization, pca, whitening, powerTrans)
%   `pca`:: 500
%       set to Inf to disable pca
%   `pcaNumSamples`:: 10^5
%       set to Inf to include all samples 
%   `whiten`:: true
%       set to false to diable whitening
%   `powerTrans`:: 0.5
%       set to 1 to disable power transform

if nargin<2 || isempty(model),
    model = 'imagenet-vgg-m';
end
if nargin<1 || isempty(imdbName),
    imdbName = 'ModelNet40v1';
end
if ischar(model), 
    modelName = model; 
    net = [];
else
    modelName = 'NoName';
    net = model;
end

% default options
opts.aug = 'none';
opts.gpuMode = false;
opts.numWorkers = 12;
opts.restart = false;
opts.readOp = @imread_255;
opts.normalization = true;
opts.pca = 500;
opts.pcaNumSamples = 10^5;
opts.whiten = true;
opts.powerTrans = 0.5;
[opts,varargin] = vl_argparse(opts,varargin);

% -------------------------------------------------------------------------
%                                                                 Get imdb
% -------------------------------------------------------------------------
imdb = get_imdb(imdbName);
[imdb.images.id,I] = sort(imdb.images.id);
I = I(1:11808);%first 10 classes
imdb.images.name = imdb.images.name(I);
imdb.images.class = imdb.images.class(I);
imdb.images.set = imdb.images.set(I);
if isfield(imdb.images,'sid'), imdb.images.sid = imdb.images.sid(I); end
nImgs = numel(imdb.images.name);

% -------------------------------------------------------------------------
%                                                                CNN Model
% -------------------------------------------------------------------------
if isempty(net),
    netFilePath = fullfile('/media/DATA/mvcnn','models', [modelName '.mat']);
    % download model if not found
    if ~exist(netFilePath,'file'),
        fprintf('Downloading model (%s) ...', modelName) ;
        vl_xmkdir(fullfile('/media/DATA/mvcnn','models')) ;
        urlwrite(fullfile('http://maxwell.cs.umass.edu/deep-shape-data/models', ...
            [modelName '.mat']), netFilePath) ;
        fprintf(' done!\n');
    end
    net = load(netFilePath);
    net = convert_net_format(net,'old'); 
end

% use gpu if requested and possible
if opts.gpuMode,
    if gpuDeviceCount()==0,
        fprintf('No supported gpu detected! ');
        reply = input('Continue w/ cpu mode? Y/N [Y]:','s');
        if ~isempty(reply) && reply~='Y',
            return;
        end
        opts.gpuMode = false;
    else
        gpuDevice(1);
        net = vl_simplenn_move(net,'gpu');
    end
end

% see if it's a multivew net
viewpoolIdx = find(cellfun(@(x)strcmp(x.name, 'viewpool'),net.layers));
if ~isempty(viewpoolIdx), 
    if numel(viewpoolIdx)>1, 
        error('More than one viewpool layers found!'); 
    end
    nViews = net.layers{viewpoolIdx}.stride;
    [imdb.images.sid,I] = sort(imdb.images.sid);
    imdb.images.name = imdb.images.name(I);
    imdb.images.class = imdb.images.class(I);
    imdb.images.set = imdb.images.set(I);
    imdb.images.id = imdb.images.id(I);
else
    nViews = 1;
end
nImgs = nImgs / nViews;

% response dimensions
fprintf('Testing model (%s) ...', modelName) ;
if isfield(net.layers{1},'weights'), 
  nChannels = size(net.layers{1}.weights{1},3); 
else
  nChannels = size(net.layers{1}.filters,3);  % old format
end
im0 = zeros(net.normalization.imageSize(1), ...
    net.normalization.imageSize(2), nChannels, nViews, 'single') * 255; 
if opts.gpuMode, im0 = gpuArray(im0); end
res = vl_simplenn(net,im0);
layers.name = {};
layers.sizes = [];
layers.index = [];
for i = 1:length(net.layers), 
    ires = i+1;
    [sz1, sz2, sz3, sz4] = size(res(ires).x);
%     if sz1==1 && sz2==1 && sz4==1 && isfield(net.layers{i},'name'),
%         layers.name{end+1} = net.layers{i}.name;
%         layers.index(end+1) = ires;
%         layers.sizes(:,end+1) = [sz1;sz2;sz3];
%     end
    if sum(i == 13:15) ~= 0,
        layers.name{end+1} = net.layers{i}.name;
        layers.index(end+1) = ires;
        layers.sizes(:,end+1) = [sz1;sz2;sz3];
    end
end
fprintf(' done!\n');

% -------------------------------------------------------------------------
%                                                   Load data if available
% -------------------------------------------------------------------------
% saving directory
saveDir = fullfile('/media/DATA/mvcnn','features',sprintf('%s-%s-%s', ...
    imdbName, modelName, opts.aug));
if ~opts.normalization, 
    expSuffix = 'NORM0';
else
    expSuffix = sprintf('NORM%d-PCA%d', opts.normalization, opts.pca);
end
if opts.restart,
    rmdir(saveDir,'s');
end
if nViews>1, 
    cacheDir = fullfile(saveDir,'cache','sid');
else
    cacheDir = fullfile(saveDir,'cache','id');
end
vl_xmkdir(cacheDir);
vl_xmkdir(fullfile(saveDir,expSuffix));

featCell = cell(1,numel(layers.name));
flag_found = true;
fprintf('Loading pre-computed features ... ');
for fi = 1:numel(layers.name),
    featPath = fullfile(saveDir,expSuffix,[layers.name{fi} '.mat']);
    if ~exist(featPath, 'file'), 
        flag_found = false;
        break;
    end
    fprintf('%s ... ', layers.name{fi});
    featCell{fi} = load(featPath);
end
if flag_found, 
    fprintf('all found! \n');
    feats = struct();
    for fi = 1:numel(layers.name),
        feats.(layers.name{fi}) = featCell{fi};
    end
    return;
else
    fprintf('all/some feature missing! \n');
    clear featCell;
end

% -------------------------------------------------------------------------
%                                                    Get raw CNN responses
% -------------------------------------------------------------------------
% data augmentation
subWins = get_augmentation_matrix(opts.aug);
nSubWins = size(subWins,2);

if opts.numWorkers<=1 || opts.gpuMode, 
    poolSize = 0;
else
    poolObj = gcp('nocreate');
    if isempty(poolObj) || poolObj.NumWorkers<opts.numWorkers, 
        if ~isempty(poolObj), delete(poolObj); end
        poolObj = parpool(opts.numWorkers);
    end
    poolSize = poolObj.NumWorkers;
end
parfor_progress(nImgs); 
parfor (i=1:nImgs, poolSize)
%  for  i=1:nImgs, % if no parallel computing toolbox
    if exist(fullfile(cacheDir, [num2str(i) '.mat']),'file'),
        continue;
    end
    im = cell(1,nViews);
    for v = 1:nViews, 
        im{v} = opts.readOp(...
            fullfile(imdb.imageDir, imdb.images.name{(i-1)*nViews+v}), ...
            nChannels);
        if isfield(imdb.meta,'invert') && imdb.meta.invert,
            im{v} = 255 - im{v};
        end
    end
    
    feat = get_cnn_activations( im, net, subWins, layers, ...
        'gpuMode', opts.gpuMode);
    parsave(fullfile(cacheDir, [num2str(i) '.mat']),feat);
    
    % fprintf(' %s\n',fullfile(imdb.imageDir,imdb.images.name{(i-1)*nViews+1}));
    parfor_progress();
end
parfor_progress(0);

% -------------------------------------------------------------------------
%                               Construct feature descriptors and encoders
% -------------------------------------------------------------------------
featCell = cell(1,numel(layers.name));
for fi=1:numel(layers.name),
    featCell{fi}.x = single(zeros(nImgs*nSubWins,layers.sizes(1,fi),layers.sizes(2,fi),layers.sizes(3,fi)));
    featCell{fi}.imdb = imdb;
    featCell{fi}.modelName = modelName;
    featCell{fi}.layerName = layers.name{fi};
    if nViews>1, 
        featCell{fi}.sid = reshape(repmat([1:nImgs],[nSubWins,1]),...
            [nImgs*nSubWins,1]);
    else
        featCell{fi}.id = reshape(repmat([1:nImgs],[nSubWins,1]),...
            [nImgs*nSubWins,1]);
    end
end
fprintf('Loading raw features: \n');
for i=1:nImgs,
    if mod(i,10)==0, fprintf('.'); end
    if mod(i,500)==0, fprintf(' %4d/%4d\n', i,nImgs); end
    feat = load(fullfile(cacheDir, [num2str(i) '.mat']));
    for fi = 1:numel(layers.name),
        featCell{fi}.x((i-1)*nSubWins+(1:nSubWins),:,:,:) = ...
            squeeze(feat.(layers.name{fi}));
%             squeeze(feat.(layers.name{fi}))';
    end
end
fprintf(' %4d/%4d done! \n', nImgs,nImgs);

if opts.normalization,
    fprintf('Normalization: \n');
    nSamples = opts.pcaNumSamples;
    if nSamples >= nImgs*nSubWins,
        nSamples = nImgs*nSubWins;
    end
    for fi = 1:numel(layers.name),
        featOldDim = size(featCell{fi}.x,2);
        % TODO deal with singlularity
        if featOldDim>nSamples,
            error('No enough samples');
        end
        featNewDim = opts.pca;
        if featNewDim>=featOldDim,
            featNewDim = featOldDim;
        end
        fprintf('[%d/%d] (%s): ', fi, numel(layers.name), layers.name{fi});
        fprintf('l2 ... ');
        featCell{fi}.x = bsxfun(@rdivide,featCell{fi}.x,...
            sqrt(sum(featCell{fi}.x.^2,2)));
        fprintf('pca (%d-D=>%d-D,N=%d,whiten=%d) ... ', featOldDim, ...
            featNewDim, nSamples, opts.whiten);
        roIdx = randperm(nImgs*nSubWins);
        x = featCell{fi}.x(roIdx(1:nSamples),:);
        featCell{fi}.pcaMean = mean(x);
        [coeff,~,latent] = pca(x);
        if opts.whiten,
            coeff = bsxfun(@rdivide,coeff,sqrt(latent)');
        end
        featCell{fi}.pcaCoeff = coeff(:,1:featNewDim);
        featCell{fi}.x = bsxfun(@minus,featCell{fi}.x,featCell{fi}.pcaMean) ...
            * featCell{fi}.pcaCoeff;
        fprintf('l2 ... ');
        featCell{fi}.x = bsxfun(@rdivide,featCell{fi}.x,...
            sqrt(sum(featCell{fi}.x.^2,2)));
        featCell{fi}.powerTrans = opts.powerTrans;
        if featCell{fi}.powerTrans~=1,
            fprintf('pow (%.1f) ... ',featCell{fi}.powerTrans);
            featCell{fi}.x = (abs(featCell{fi}.x).^featCell{fi}.powerTrans)...
                .* sign(featCell{fi}.x);
        end
        fprintf('done! \n');
    end
end

% write to disk
feats = struct();
fprintf('Saving feature descriptors & encoders: ');
for fi = 1:numel(layers.name),
    fprintf('%s ... ',layers.name{fi});
    feat = featCell{fi};
    feats.(layers.name{fi}) = feat;
    save(fullfile(saveDir, expSuffix, [layers.name{fi} '.mat']), ...
        '-struct', 'feat', '-v7.3');
end
fprintf('done! \n');

end
