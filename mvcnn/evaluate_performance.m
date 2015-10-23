function evaluate_performance( imdbName, model, varargin )
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
opts.logPath = fullfile('log','eval_new.txt');
[opts,varargin] = vl_argparse(opts,varargin);

% -------------------------------------------------------------------------
%                                                                 Get imdb
% -------------------------------------------------------------------------
if ~exist('imdb','var'), 
    load('data/fc6.mat','imdb')
end
% imdb = get_imdb(imdbName);
[imdb.images.id,I] = sort(imdb.images.id);
imdb.images.name = imdb.images.name(I);
imdb.images.class = imdb.images.class(I);
imdb.images.set = imdb.images.set(I);
if isfield(imdb.images,'sid'), imdb.images.sid = imdb.images.sid(I); end
nImgs = numel(imdb.images.name);

% -------------------------------------------------------------------------
%                                                                CNN Model
% -------------------------------------------------------------------------
if isempty(net),
    netFilePath = fullfile('./data','models', [modelName '.mat']);
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
% if isfield(net.layers{1},'weights'), 
%   nChannels = size(net.layers{1}.weights{1},3); 
% else
%   nChannels = size(net.layers{1}.filters,3);  % old format
% end
% im0 = zeros(net.normalization.imageSize(1), ...
%     net.normalization.imageSize(2), nChannels, nViews, 'single') * 255; 
% if opts.gpuMode, im0 = gpuArray(im0); end

if ~exist('train_data','var'), 
    load('data/train_data.mat')
    load('data/test_data.mat')
    load('data/train_labels.mat')
    load('data/test_labels.mat')
end

im_train(1,1,:,:) = single(train_data');
res_train = vl_simplenn(net,im_train);
train_softmax_result = res_train(end).x(1,1,:,:);
[~,train_result] = max(train_softmax_result);
accuTrain = sum(train_result(:) == train_labels)/numel(train_labels)


im_test(1,1,:,:) = single(test_data');
res_test = vl_simplenn(net,im_test);
test_softmax_result = res_test(end).x(1,1,:,:);
[~,test_result] = max(test_softmax_result);
accuTest = sum(test_result(:) == test_labels)/numel(test_labels)

fid = fopen(opts.logPath,'a+');
fprintf(fid, '(%s) -- Classification\n', datestr(now));
fprintf(fid, '\tmodel: %s\n',modelName);
fprintf(fid, '\taccuracy (train): %g%%\n',accuTrain*100);
fprintf(fid, '\taccuracy (test): %g%%\n',accuTest*100);
fclose(fid);