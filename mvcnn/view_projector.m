addpath(genpath('./'))
if ~exist('MAP','var')
    load('/media/DATA/mvcnn/data/LEARN_MODEL_DIM-class_split-test.mat')
end

%% Generate set of poses around unit circle

ncent = 12;% mapping centers, Need to know if this number affects the performance???????

icent=((1:ncent)-1)'*2*pi/ncent;

cents = [cos(icent) sin(icent)];

n_initial_poses = 12; 

initial_poses=((1:n_initial_poses)-1)'*2*pi/n_initial_poses;

angle_from_index = @(indx)((indx-1)*2*pi/n_initial_poses);

x = [cos(initial_poses) sin(initial_poses)];

G = full(direct_gaussian_RBF(x,cents));%G nX8

 

%% use cose tensor CORE_TENSOR and test image to generate set of

%coefficients

%this order verified from getCoeffTensor_cell function

%CF is wide matrix before being vectorized, as follows

% CORE_TENSOR = MAP.S*MAP.V;
% 
% Basis = MAP.U;

D = 1000;

CORE_TENSOR = MAP.S(1:D,1:D)*MAP.V(1:D,:);
Basis = MAP.U(:,1:D);

clear MAP;

 CORE_TENSOR = reshape(CORE_TENSOR' ,[],size(G,2),size(CORE_TENSOR,1));

 CORE_TENSOR_pose = tmul(CORE_TENSOR,G',2);

 % % multiply the initial pose vectores G by the core tensoe (A\epsi(x))

clear CORE_TENSOR

CORE_TENSOR_pose = mat2cell(CORE_TENSOR_pose,size(CORE_TENSOR_pose,1),ones(n_initial_poses,1),size(CORE_TENSOR_pose,3));

CORE_TENSOR_pose = cellfun(@(p)(squeeze(p)'),CORE_TENSOR_pose,'UniformOutput',false);