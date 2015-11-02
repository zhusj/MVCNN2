function [class_recog_rate,pose_recog_rate] = pose_category_estimation(pd,run_stamp,p_trained_data_file,p_test_data_file,dim,p_classifier,p_save_file)

%% adjust process parameters
load('sampling_params');
if ~exist('P_CLASS_REFINEMENT','var')
    P_CLASS_REFINEMENT = 0;
end
if ~exist('P_GRADIENT_OPTMZ','var')
    P_GRADIENT_OPTMZ = 0;
end

if ~exist('P_N_SMPLNG_ITRTNs','var')
    P_N_SMPLNG_ITRTNs = 200;
end

if ~exist('P_PARLL_TEST','var')
    P_PARLL_TEST = 5;
end

if ~exist('P_PLOT_SAMPLES','var')
    P_PLOT_SAMPLES = 0;
end

if ~exist('P_POSE_REFINEMENT','var')
    P_POSE_REFINEMENT = 0;
end

if ~exist('P_TEST_RATE','var')
    P_TEST_RATE = 1;
end
if ~exist('P_USE_DRUR_REG','var')
    P_USE_DRUR_REG = 0;
end

if ~exist('P_START_INDEX','var')
    P_START_INDEX = 1;
end
if ~exist('P_SAMPLING_ON','var')
    P_SAMPLING_ON=1;
end

if ~exist('P_CATEGORICAL_POSE','var')
    P_CATEGORICAL_POSE=1;
end
if ~exist('P_GRAD_COST_FUN','var')
    P_GRAD_COST_FUN=2;
end
if ~exist('P_GRAD_N_INIT_PTS','var')
    P_GRAD_N_INIT_PTS=30;
end


if P_USE_DRUR_REG
    classifier = 'SVM';
else
    classifier = 'RfC';%regression for classification
end
if ~isempty(strfind(p_classifier,'SVM'))  % Using SVM Classifier
    classifier = 'SVM';
end
if ~isempty(strfind(p_classifier,'KNN'))  % Using SVM Classifier
    classifier = 'KNN';
end
%%
load(p_trained_data_file);
MAP.classification = 1;
MAP = MAP.setDimension(pd);
T_train = MAP.T;

%T_train = T_train(:,2:end);



%% SET T and learn RKHS function from T to CF 
if ~P_USE_DRUR_REG
    % Direct
    if length(unique(train_labels))==length(train_labels)
        T_cents = T_train;
    else
        T_cents = [];
        for i = unique(train_labels')
            indx = (train_labels==i);
            t = T_train(indx,:);
            T_cents = [T_cents;mean(t)];
        end
    end
    T_CF_map = GRN_RKHS;
    %T_CF_map  = T_CF_map.learn_regularized(T_train,MAP.Orig_X,3*numel(unique(train_labels)),-1);%,T_cents);%size(tSeqs,1),0);
    T_CF_map  = T_CF_map.learn_regularized(T_train,MAP.Orig_X,30,-0.05);%,T_cents);%size(tSeqs,1),0);
    %T_cents = T_CF_map.Cents;
else
    %Using DRUR
    drur_done =0;
    if exist(sprintf('DRUR_data_%s.mat',run_stamp),'file')
        load(sprintf('DRUR_data_%s',run_stamp),'T_train','T_cents','ff','T_CF_map');
        if size(T_train,2) == pd
            drur_done = 1;
        end
    end

    if ~drur_done
        [XX,ff,FF,E] = drur(MAP.Orig_X,T_train,{'rbf',15,[],0.1},{'rbf',20,[],0.1},[],20);
        T_train = XX{end}; ff = ff{end};FF = FF{end};

        T_cents = [];
        for i = unique(train_labels')
            indx = (train_labels==i);
            t = T_train(indx,:);
            T_cents = [T_cents;mean(t)];
        end

        % %load('DRUR_data_indeces-train');
        % %Tdrur_test = rbf(A_test,FF);
        T_CF_map.predict_regularized = @(p_t) rbf(p_t,ff);
        %DRUR_CF_T_map.predict_regularized = @(p_cf) rbf(p_cf,FF);
        save(sprintf('DRUR_data_%s',run_stamp),'T_train','T_cents','ff','FF','T_CF_map');
    end
end

T_stat.dim_max = max(T_train);
T_stat.dim_min = min(T_train);
T_stat.dim_mean = mean(T_train);
T_stat.dim_var = var(T_train);

%clear T; clear T_cents;
%% learn SVM model for latent space
SVM_model = create_SVM_model(T_train, train_labels,1);
%%if strcmp(classifier,'SVM')
% if exist(sprintf('SVM_model%s.mat',run_stamp),'file')
%      load(sprintf('SVM_model%s',run_stamp));
% else
%     SVM_model = create_SVM_model(T_train, train_labels,1);
%     %SVM_model = tg_SVMtrain(T_train, train_labels, 8);
%     %SVM_model = ovrtrain(train_labels, T_train,'-t 0 -q');%linear SVM model
%     %SVM_model = svmtrain(T_train, train_labels,0);
%      save(sprintf('SVM_model%s',run_stamp),'SVM_model');
%  end
%%end

%% LOAD Test images and its labels
load(p_test_data_file);
N = size(ACTUAL_IMAGES,1);

if exist('P_TEST_RATE','var') && P_TEST_RATE>1
    i = 1:P_TEST_RATE:N;
    %i = randsample(N,ceil(N/P_TEST_RATE),false);
    l = false(N,1);
    l(i) = true;
    ACTUAL_IMAGES = ACTUAL_IMAGES(l,:);    
else
    l = true(N,1);
end

ACTUAL_CLASS_LABELS = ACTUAL_CLASS_LABELS(l,:);
ACTUAL_INSTANCE_LABELS = ACTUAL_INSTANCE_LABELS(l,:);
ACTUAL_HIGHT_LABELS = ACTUAL_HIGHT_LABELS(l,:);
if exist('ACTUAL_SCALE_LABELS','var') %for the RGBD dataset doesn't has scale
    ACTUAL_SCALE_LABELS = ACTUAL_SCALE_LABELS(l,:);
else
    ACTUAL_SCALE_LABELS = -1*ones(length(ACTUAL_HIGHT_LABELS),1) ;
end
ACTUAL_POSE_LABELS = ACTUAL_POSE_LABELS(l,:);

N = size(ACTUAL_IMAGES,1);
    
%% SET UP POSE IMAGE MAP
if P_CATEGORICAL_POSE
	nposes = max(ACTUAL_POSE_LABELS);
	iposes = ((1:nposes)-1)'*2*pi/nposes;
    pose_classify = @(p_poses)knnclassify(mod(p_poses,2*pi),[iposes;2*pi] ,[1:nposes 1]');
else
    pose_classify = @(p_poses)mod(p_poses,2*pi)*180/pi;
end


load('all_GCF','ncents','cents');
if ~exist('ncents','var')
    ncents = 8;%size(tImgs,5);
end
icents=((1:ncents)-1)'*2*pi/ncents;
if ~exist('cents','var')
    cents = [cos(icents) sin(icents)]; 
end

poses_stat.dim_max = max(icents);
poses_stat.dim_min = min(icents);
poses_stat.dim_mean = mean(icents);
poses_stat.dim_var = max(icents);

POSE_img_map = GRN_RKHS;
POSE_img_map = POSE_img_map.set_cents(cents);

%% START SAMPLING PROCESS

init_prtcl_stats.dim_max = [T_stat.dim_max poses_stat.dim_max];
init_prtcl_stats.dim_min = [T_stat.dim_min poses_stat.dim_min];
init_prtcl_stats.dim_mean = [T_stat.dim_mean poses_stat.dim_mean];
init_prtcl_stats.dim_var = [T_stat.dim_var poses_stat.dim_var];

n_parll_imgs = min(P_PARLL_TEST,N);%NUmber of images processed in parallel
total_results = zeros(N,8);%cell(size(ACTUAL_IMAGES,1)/n_parll_imgs,1);
sampling_time = zeros(N,1);
best_prtcls_all = zeros(N,pd+1);

if size(T_cents,1)<=20
    init_T = T_cents;
else
    [~,init_T] = kmeans(MAP.T,20);
end
%init_T = T_cents;%mean(MAP.T);

size(init_T)
t_dimension = size(init_T,2);
%T = T(randi(size(T,1),20,1),:);

n_initial_poses = 8;
initial_poses=((1:n_initial_poses)-1)'*2*pi/n_initial_poses;
pose_dim = size(initial_poses,2);
%poses = [cos(icent) sin(icent)]; 
%test_ind = randi(size(ACTUAL_IMAGES,1),n_parll_imgs,1);
%matlabpool;
if exist('P_GRADIENT_OPTMZ','var') && P_GRADIENT_OPTMZ
    R = reshape(T_CF_map.Map',[],POSE_img_map.n_cents,T_CF_map.n_cents);
    ephi = @(t)[{T_CF_map.get_kernel(t)'} {directs_grbf_jacobian(t,T_CF_map.Cents)}];
    epsi = @(x)[{POSE_img_map.get_kernel(x)'} {directs_grbf_jacobian(x,POSE_img_map.Cents)}];
    
    CostFun = P_GRAD_COST_FUN; % 0 (Euc), 1 (NCC without subtracting mean), 2 (NCC with subtracting mean)
    minMethod = 2;  %  minMethod =1 or 2 or 3
    %t_dimension = size(best_T,2);
    %optns = optOptions;
    %optns([1,2,3]) = [1,1E-6,1E-6];
    
    
else
    P_GRADIENT_OPTMZ=0;
end
bunch_size = 200;%length(ACTUAL_CLASS_LABELS)/20;
indx = P_START_INDEX;
while indx <= N
    test_ind = indx:min(indx+n_parll_imgs-1,N)
    goals  = ACTUAL_IMAGES(test_ind,:);
    %goals = cell2mat(goals')';
    
    results = zeros(length(test_ind),8);
    results(:,1) = ACTUAL_CLASS_LABELS(test_ind,:);
    results(:,2) = ACTUAL_INSTANCE_LABELS(test_ind,:);
    results(:,3) = ACTUAL_HIGHT_LABELS(test_ind,:);
    results(:,4) = ACTUAL_SCALE_LABELS(test_ind,:);
    results(:,5) = ACTUAL_POSE_LABELS(test_ind,:);
    
    plot_handles = zeros(length(test_ind));
    %
    %fn = @(x,y)dist2(MAP.predict_regularized([cos(x) sin(x)]),y)
    % fn = @(t,bases)dist2(POSE_img_map.predict_regularized([cos(p_ind), sin(p_ind)]),bases)
    % EST_SAMPLE = ParticleFilter(fn);

    % create particles from combining initial Ts and initial poses
    %t_dimension = size(init_T,2);
    t_ind = 1:size(init_T,1);
    p_ind = 1:size(initial_poses,1);
    %p_ind = round(mean(1:size(initial_poses,1)));
    goal_ind = 1:size(goals,1);
    [mex_t,mex_p,mex_g] = meshgrid(t_ind,p_ind,goal_ind);
    initial_prtcles = [init_T(mex_t(:),:),initial_poses(mex_p(:),:)];
    initial_prtcl_goal_map = mex_g(:);
    
    % set up the energy function
    fn = @(p_prtcls,p_prtcl_goal_map)energy_fun(goals,p_prtcls,p_prtcl_goal_map,T_CF_map,POSE_img_map,pose_dim);
    %% two way sampling
    if P_SAMPLING_ON
        if P_PLOT_SAMPLES 
            for i = 1:length(test_ind)
                plot_handles(i) = plotting(T_train,train_labels);
                title(sprintf('The goal class is %d',results(i,1)));
            end
        end
        % start particle filter for both pose and T
        EST_SAMPLE = ParticleFilter(fn,P_N_SMPLNG_ITRTNs,init_prtcl_stats,P_PLOT_SAMPLES,plot_handles,{1:t_dimension});
        tic
        [best_prtcls,err,EST_SAMPLE]=EST_SAMPLE.optimize_function(initial_prtcles,initial_prtcl_goal_map);%?????????????????
        best_prtcl_goal_map = 1:size(best_prtcls,1);
        sampling_time(test_ind) = toc;

        best_T = best_prtcls(:,1:t_dimension);
        
        if strcmp(classifier,'SVM')
            est_labels = ovrpredict(results(:,dim),best_T,SVM_model);
            %est_labels  = tg_SVMtest(best_T, results(:,dim), SVM_model);
        else
            if strcmp(classifier,'RfC') %the default
                est_labels = MAP.Predict(best_T,'latent');
                [~,est_labels] = max(est_labels ,[],2);
            else
                est_labels = knnclassify(best_T ,T_train,train_labels);
            end
        end
        results(:,end-2) = est_labels;

        best_poses = best_prtcls(:,end-pose_dim+1:end);
        best_poses = mod(best_poses,2*pi);%(n_initial_poses)*2*pi/n_initial_poses);
        results(:,end-1) = pose_classify(best_poses);
        if results(:,end-1)~=results(:,5)
            disp 'error in pose';
        end
    end
    
     %% refine the pose estimation
     if  P_POSE_REFINEMENT
        %p_ind = 1:size(initial_poses,1);
        goal_ind = 1:size(goals,1);
        %[mex_p,mex_g] = meshgrid(p_ind,goal_ind);
        %initial_prtcles = [initial_poses(mex_p(:),:)];
        %initial_prtcl_goal_map = mex_g(:);
        initial_prtcles  = best_prtcls(:,end);
        initial_prtcl_goal_map  = goal_ind;
        
        CF = squeeze(T_CF_map.predict_regularized(best_T));
        if size(T_train,1) >1
            CF = reshape(CF',[],POSE_img_map.n_cents,size(best_T,1));
        else
            CF = reshape(CF,[],POSE_img_map.n_cents);
        end
        
        dist_poles =  energy_fun_pose(goals,CF,best_poses,initial_prtcl_goal_map,POSE_img_map);
        dist_poles = [dist_poles,energy_fun_pose(goals,CF,best_poses+pi,initial_prtcl_goal_map,POSE_img_map)];
        
        [~,best_poles]= min(dist_poles,[],2);
        if sum(best_poles)>numel(best_poles)
            disp 'pose_refined'
        end
        best_poses = best_poses + (best_poles-1)*pi;
%         fn = @(p_prtcls,p_prtcl_goal_map)energy_fun_pose(goals,CF,p_prtcls,p_prtcl_goal_map,POSE_img_map,1);
%         EST_SAMPLE = ParticleFilter(fn,100,poses_stat);
%         [best_poses,err]=EST_SAMPLE.optimize_function(initial_prtcles,initial_prtcl_goal_map);%?????????????????
        
         
         best_poses = mod(best_poses,2*pi);%(n_initial_poses-0.5)*2*pi/n_initial_poses);
         best_prtcls(:,end) = best_poses;
         results(:,end-1) = pose_classify(best_poses);
     end
    

    %% refine category estimation
    if P_CLASS_REFINEMENT
        fn = @(p_prtcls,p_prtcl_goal_map)energy_fun_class(goals,best_poses,p_prtcls,p_prtcl_goal_map,T_CF_map,POSE_img_map);
        EST_SAMPLE = ParticleFilter(fn,150,T_stat,P_PLOT_SAMPLES,plot_handles,{[1:t_dimension]});
        [best_T,err]=EST_SAMPLE.optimize_function(best_T,(1:size(goals,1))');%?????????????????
        
        best_prtcls(1,end-1) = best_T;
        if strcmp(classifier,'SVM')
            est_labels = ovrpredict(results(:,1),best_T,SVM_model );
        else
            if strcmp(classifier,'RfC') %the default
                est_labels = MAP.Predict(best_T,'latent');
                [~,est_labels] = max(est_labels ,[],2);
            else
                est_labels = knnclassify(best_T,T_train,train_labels);
            end
        end
        results(:,end-2) = est_labels;
    end
    %% Gradient optimization
    if P_GRADIENT_OPTMZ
        if ~exist('best_prtcls','var')
            best_prtcls = initial_prtcles;
            best_prtcl_goal_map = initial_prtcl_goal_map;
        end
        best_T = zeros(size(goals,1),t_dimension);
        best_poses = zeros(size(goals,1),1);
        err = [];%zeros(size(goals,1),1);
        for i = 1:size(goals,1)
            I = goals(i,:)';
            gradient_optimize_fn = @(y) PredictTandX( y',I, R,t_dimension,minMethod, ephi,epsi,CostFun)';
            %fun = @(x) optfunc(x',I,R,t_dimension,ephi,epsi,CostFun,1)';
            %grad = @(x) gradfunc(x',I,R,t_dimension,ephi,epsi,CostFun,1)';
            
            ind = (best_prtcl_goal_map==i);
            start_prtcls = best_prtcls(ind,:);
            start_prtcls_goal_map = best_prtcl_goal_map(ind);
            
            start_prtcls = mat2cell(start_prtcls,ones(length(start_prtcls_goal_map),1));%,size(start_prtcls,2));
            start_prtcls_goal_map = num2cell(start_prtcls_goal_map);
            out = cellfun(fn,start_prtcls,start_prtcls_goal_map );
            [~,ind] = sort(out);
            ind = ind(1:min(P_GRAD_N_INIT_PTS,length(ind)));
            %ind = unique(randsample(length(start_prtcls_goal_map),20,true,exp(-5*out)));
            start_prtcls = start_prtcls(ind);
            %%%%start_prtcls_goal_map = start_prtcls_goal_map(ind);

%             [xo,~,~,fvalo]=scg(fun,[ti,thetai],optns,grad);
%             to = xo(1:m); thetao = xo(m+1:end);
            out = cellfun(gradient_optimize_fn, start_prtcls, 'UniformOutput',false); %start optimization
            out = cell2mat(out);
            [~,ind] = min(out(:,end));
            
            best_T(i,:) = out(ind,1:t_dimension);
            best_poses(i) = out(ind,t_dimension+1:t_dimension+1);
            err =[err;out(ind,end)];
        end
        best_prtcls = [best_T,best_poses];
        
        %best_prtcls(1,end-1) = best_T;
        if strcmp(classifier,'SVM')
            est_labels = ovrpredict(results(:,1),best_T,SVM_model );
        else
            if strcmp(classifier,'RfC') %the default
                est_labels = MAP.Predict(best_T,'latent');
                [~,est_labels] = max(est_labels ,[],2);
            else
                est_labels = knnclassify(best_T,T_train,train_labels);
            end
        end
        
        best_poses = mod(best_poses,2*pi);%(n_initial_poses-0.5)*2*pi/n_initial_poses);
        best_prtcls(:,end) = best_poses;
        
        %nnz1 = nnz(results(:,[dim 5])-results(:,[6 7]));
        
        %disp 'results before gradient'
        %results
        results(:,end-1) = pose_classify(best_poses);
        results(:,end-2) = est_labels;
        %disp 'results after gradient'
        %results
%         nnz2 = nnz(results(:,[dim 5])-results(:,[6 7]));
%         
%         if nnz2 < nnz1
%             disp 'gradient improved the result'
%         end
    end
%%
    
    %results(:,end-1) = knnclassify(poses,icents,(1:ncents)');

    results(:,end) = err;
    results(:,[dim,end-2,end-3,end-1,end])
    %total_results{ceil(indx/n_parll_imgs)} = results;
    total_results(test_ind,:) = results;
    best_prtcls_all(test_ind,:) = best_prtcls;
    clear best_prtcls
    clear results
    
    diff = abs(total_results(1:max(test_ind),[dim 5])-total_results(1:max(test_ind),[6 7]));
    class_recog_rate = 1- (nnz(diff(:,1))/size(diff,1))
    if P_CATEGORICAL_POSE
        pose_recog_rate = 1- (nnz(diff(:,2))/size(diff,1))
    else
        pose_recog_rate = 1-(sum(min(diff(:,2),360-diff(:,2))/180)/size(diff,1))
    end

    
    if any(mod(test_ind,bunch_size)==0)
        disp 'here is new bunch';
        save(sprintf('temp_results_%s',run_stamp),'indx','total_results','best_prtcls_all','class_recog_rate', 'pose_recog_rate');
        disp 'temp_saved'
    end
    indx = indx + n_parll_imgs;
    
end
%%
%results = total_results;%cell2mat(total_results');
fileName = sprintf('%s_%s_%s',mfilename, run_stamp);
save(fileName ,'total_results','class_recog_rate', 'pose_recog_rate','classifier','pose_classify','best_prtcls_all','pose_dim','P_*');
%pose_dim = 1;
best_T = best_prtcls_all(:,1:end-pose_dim);
n_t = size(best_T,1);

SVM.class_estimate = ovrpredict(ones(n_t,1),best_T,SVM_model);
diff = total_results(:,dim) - SVM.class_estimate ;
SVM.class_recog_rate = 1- (nnz(diff)/n_t);

KNN.class_estimate = knnclassify(best_T,T_train,train_labels);
KNN.class_recog_rate  = 1- (nnz(total_results(:,dim) - KNN.class_estimate)/n_t) ;

save(fileName ,'SVM','KNN','-append');
%if strcmp(classifier,'RfC')
    est_labels = MAP.Predict(best_T,'latent');
    [~,est_labels] = max(est_labels ,[],2);
    RfC.class_estimate = est_labels;
    diff = total_results(:,dim)-RfC.class_estimate;
    RfC.class_recog_rate = 1- (nnz(diff)/n_t);
    
% else
%     RfC.class_estimate = [];
%     RfC.class_recog_rate = -1;
% end
save(fileName ,'RfC','-append');

classes = unique(total_results(:,dim))';
class_confusion = zeros(max(classes),max(classes));
for i = classes
    ind = (total_results(:,1)==i);
    for j = classes
        class_confusion(i,j) = sum(total_results(ind,end-2)==j)/sum(ind);
    end
end
h1=figure; imagesc(class_confusion);
print(h1, '-dpng', sprintf('class_confusion_%s',run_stamp));
save(fileName ,'class_confusion','-append');

%% pose numbers
est_poses = best_prtcls_all(:,end);
est_poses = mod(est_poses,2*pi);
if P_CATEGORICAL_POSE
    actual_poses = iposes(ACTUAL_POSE_LABELS);
    diff_poses = abs(est_poses - actual_poses);
    diff_poses2 = diff_poses;
    diff_poses2(ACTUAL_POSE_LABELS==1) = min(diff_poses(ACTUAL_POSE_LABELS==1),abs(est_poses(ACTUAL_POSE_LABELS==1)-2*pi));
else
    actual_poses = ACTUAL_POSE_LABELS;
    diff_poses = abs(est_poses - actual_poses);
    diff_poses2 = min(diff_poses,360-diff_poses);
end

diff_less_45 = (diff_poses2<(pi/4));
pose_rate_less_45 = sum(diff_less_45)/length(diff_less_45)
diff_less_225 = (diff_poses2<(pi/8));
pose_rate_less_225 = sum(diff_less_225)/length(diff_less_225)


save(fileName ,'pose_rate_less_225','pose_rate_less_45','-append');


% save data to file
%if p_save_file
    %save(sprintf('%s_%s',mfilename, run_stamp),'total_results','class_recog_rate', 'pose_recog_rate','classifier','best_prtcls_all','SVM','RfC','pose_rate_less_225','pose_rate_less_45','icents');
%end

if P_CATEGORICAL_POSE
    %nposes = length(iposes);
    pose_confusion = zeros(nposes,nposes);
    for i = 1:nposes
        ind = (total_results(:,end-3)==i);
        for j = 1:nposes
            pose_confusion(i,j) = sum(total_results(ind,end-1)==j)/sum(ind);
        end
    end
    h=figure; imagesc(pose_confusion);
    print(h, '-dpng', sprintf('pose_confusion_%s',run_stamp));

    if ~mod(nposes,2)
        antipodal_pose_recog_rate=(sum(diag(pose_confusion)) +sum(diag(pose_confusion,nposes/2))+sum(diag(pose_confusion,-nposes/2)))/nposes;
    else
        antipodal_pose_recog_rate = -1;
    end
    
    save(fileName ,'pose_confusion','iposes','nposes','-append');
end

%save(sprintf('%s_%s',mfilename, run_stamp),'pose_confusion','class_confusion','antipodal_pose_recog_rate','-append');

end

%% %%%%%%%%%%%%%%%%%%%%

function parsave(fname, total_results,indx)
    save(fname, 'indx','total_results');
end


