% function run_experiments()

setup;

logPath = fullfile('log','eval_new.txt'); 
skipTrain = false; 
skipEval = false; 
ex = struct([]);

% experiment unit -- fine-tuning MVCNN using 12 views w/ upright assumption
ex(end+1).trainOpts = struct(...
    'baseModel', '', ... 
    'dataset', 'ModelNet40v1', ...
    'batchSize', 10, ... 
    'aug', 'none', ...
    'numEpochs', 10, ...
    'gpuMode', true, ...
    'multiview', false, ...
    'viewpoolLoc', 'conv5', ...
    'learningRate', [1e-4*ones(1,3) 3e-5*ones(1,3) 1e-5*ones(1,3) 1e-6*ones(1,3)]*10, ...%set1: [1e-4*ones(1,3) 3e-5*ones(1,3) 1e-5*ones(1,3) 1e-6*ones(1,3)], ...    set2: logspace(-2, -4, 60) ;set3: logspace(-3, -4, 10) ;
    'momentum', 0.5, ...
    'includeVal', true, ...
    'useUprightAssumption', false); 
ex(end).featOpts = struct(...
    'dataset', 'ModelNet40v1', ...
    'aug', 'none', ...
    'gpuMode', true, ...
    'numWorkers', 6);
ex(end).claOpts = struct(...
    'feat', 'relu7'); 
ex(end).retOpts = struct(...
    'feat','relu7',...
    'gpuMode', true, ...
    'numWorkers', 6); 


for i=1:length(ex), 
    % ---------------------------------------------------------------------
    %                                                    train / fine-tune 
    % ---------------------------------------------------------------------
    if isfield(ex(i),'trainOpts') && ~skipTrain, 
        trainOpts = ex(i).trainOpts;
%         prefix = sprintf('BS%d_AUG%s', trainOpts.batchSize, trainOpts.aug);
        prefix = 'pose_10_epochs_10_times_leariningRate_batchSize_10_W_d_500_SVM_fc7';
        if isfield(trainOpts,'multiview') && trainOpts.multiview, 
            prefix = sprintf('%s_MV%s',prefix,trainOpts.viewpoolLoc);
        end
        modelName = sprintf('%s-finetuned-%s-%s', trainOpts.baseModel, ...
            trainOpts.dataset, prefix);
%         modelName = sprintf('%s-finetuned-%s-%s', trainOpts.baseModel, ...
%             trainOpts.dataset, 'branch');
        trainOpts.prefix = prefix;
        if ~exist(fullfile('/media/DATA/mvcnn/data','models',[modelName '.mat']),'file'),
            net = run_train(trainOpts.dataset, trainOpts);
            save(fullfile('/media/DATA/mvcnn/data','models',[modelName '.mat']),'-struct','net');
        end
        if isfield(ex(i),'featOpts'), ex(i).featOpts.model = modelName; end
    end
    
    % ---------------------------------------------------------------------
    %                                                     compute features 
    % ---------------------------------------------------------------------
    clear feats;
    if isfield(ex(i),'featOpts') && ~skipEval,
        featOpts = ex(i).featOpts;
        featDir = fullfile('/media/DATA/mvcnn/data', 'features', ...
            [featOpts.dataset '-' featOpts.model '-' featOpts.aug], 'NORM0');
        if exist(fullfile(featDir, 'prob.mat'),'file'), % supposedly the last
            fprintf('Existing descriptors found at %s \n', featDir);
        end
        if ~isfield(featOpts, 'logPath'), featOpts.logPath = logPath; end
%         feats = imdb_compute_cnn_features(featOpts.dataset, featOpts.model, ...
%             'normalization', false, featOpts);
        evaluate_performance(featOpts.dataset, featOpts.model, ...
            'normalization', false, featOpts);
    end
    
    % ---------------------------------------------------------------------
    %                                            classification evaluation
    % ---------------------------------------------------------------------
%     if isfield(ex(i),'claOpts') && exist('feats','var') && ~skipEval, 
%         claOpts = ex(i).claOpts;
%         evalClaPath = fullfile(featDir,claOpts.feat,'evalCla.mat');
%         if exist(evalClaPath, 'file'), 
%             fprintf('Classification evaluated before at %s \n', evalClaPath);
%         else
%             if ~isfield(claOpts, 'log2c'), claOpts.log2c = [-8:4:4]; end
%             if ~isfield(claOpts, 'cv'), claOpts.cv = 2; end
%             if ~isfield(claOpts, 'logPath'), claOpts.logPath = logPath; end
%             evaluate_classification(feats.(claOpts.feat), ...
%                 'predPath', evalClaPath, ...
%                 claOpts);
%         end
%     end
%     
%     % ---------------------------------------------------------------------
%     %                                                 retrieval evaluation 
%     % ---------------------------------------------------------------------
%     if isfield(ex(i),'retOpts') && exist('feats','var') && ~skipEval, 
%         retOpts = ex(i).retOpts;
%         evalRetPath = fullfile(featDir,retOpts.feat,'evalRet.mat');
%         if exist(evalRetPath, 'file'), 
%             fprintf('Retrieval evaluated before at %s \n', evalRetPath);
%         else
%             if ~isfield(retOpts, 'logPath'), retOpts.logPath = logPath; end
%             [res,info] = retrieve_shapes_cnn([],feats.(retOpts.feat),retOpts);
%             save(evalRetPath,'res','info');
%         end
%     end
    
end

% delete parallel pool if there is one
delete(gcp('nocreate'));
