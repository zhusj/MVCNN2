% learn Kernel PLS mapping: implementatin to what in paper Kenrk Partial
% Least square regression in RKHS.
% input:
% D : original data matrix Nxd  N: is the number of samples, d: is the dimensionality of the input images
% X : embedding space representation. Nxe matrix N: number of samples, e: the dimesnionality of the embedding space
% cent: centers: exm e: the dimesnionality of the embedding space

classdef KPLS
    properties
        X
        mean_orig_X 
        Y
        mean_orig_Y 
        K
        classification
    end
    
    properties(SetAccess = private)
        D_
        T_
        U_
    end
    properties(Dependent = true, SetAccess= private)
%         R
        R2
%         M
        M2
    end
%     
    methods
        
        
        %% learn/predict functions
        function [r,l]=Predict(obj,Tt,p_input_type)
            if strcmp(p_input_type,'latent')
                r = Tt*obj.T'*obj.Y;
            else
                if strcmp(p_input_type,'input')
                    %Tt = bsxfun(@minus,Tt,obj.mean_orig_X);
                    %Kt = get_gram_kernel(Tt,obj.X);
                    Kt = get_gram_kernel(Tt,obj.Orig_X);
                    %Kt = At*A_train';
                    [nt,n] = size(Kt);
                    wt = ones(nt,1);
                    w = ones(n,1);
                    wt= wt*w';
                    %w = eye(n)-w/n;
                    Kt = (Kt-wt/n*obj.K)*(eye(n)-w*w'/n);
                    r = Kt*obj.M2;
                end
            end
            if obj.classification
                [~,l] = max(r,[],2);
            else
                l = r+repmat(obj.mean_orig_Y,size(r,1),1);
            end
            
        end
        
        function r=get_latent(obj,Xt)
            %Xt = bsxfun(@minus,Xt,obj.mean_orig_X);
            %Kt = get_gram_kernel(Xt,obj.X);
            Kt = get_gram_kernel(Xt,obj.Orig_X);
            %Kt = At*A_train';
            [nt,n] = size(Kt);
            wt = ones(nt,1);
            w = ones(n,1);
            wt= wt*w';
            %w = eye(n)-w/n;
            Kt = (Kt-wt/n*obj.K)*(eye(n)-w*w'/n);
            
            r = Kt*obj.R2;
        end
        %%%%%%%%%%% learning function %%%%%%%%%%%%%%%%%%
        function [MAP,Y1,normK] = Learn(MAP,X,Labels,dim,stp,p_classification)
            
            if ~exist('dim','var') || isempty(dim) || dim<1 || (dim>=min(size(X))) 
                dim = min(size(X));
            else
                dim = dim+1;
            end

            if ~exist('stp','var') || stp==0 || isempty(stp)
                stp = 1E-4;
            end
            
            if ~exist('p_classification','var')
                p_classification = 1;
            end
            
            MAP.classification = p_classification;
            %D = 101;
            n = size(X,1);

            %% normalizing/centeralizing data
            if p_classification% && (numel(Labels)==length(Labels))
                L = get_matrix_label(Labels);
            else
                L = Labels;
            end
            MAP.mean_orig_Y=mean(L);
            L = bsxfun(@minus,L, MAP.mean_orig_Y);
            MAP.Y = L;
%             sgma = sqrt(var(L));
%             L =  bsxfun(@rdivide,L,sgma);

            K = get_gram_kernel(X); %K = A*A';
            w = ones(n,1);
            w= w*w';
            w = eye(n)-w/n;
            K = w*K*w';%normalizing K using equation (13) in the paper
            MAP.K = K;
            normK = norm(K);
            MAP.mean_orig_X = mean(X);
            MAP.X = bsxfun(@minus, X, MAP.mean_orig_X);

            %% Iterating
            thresh = 1E-8;
            indx = 1;
            T = [];
            U = [];
            C = [];
            %W = [];
            S = [];
            flag=0;
            reason = 1; % indicate the reason of stopping the iteration.
            while indx<dim
                t_old = zeros(n,1);
                t = ones(n,1);
                u = rand(n,1);
                u_old = -1*u;

                itrns = 0;
                while norm(abs(t)-abs(t_old))>thresh %|| norm(u-u_old)>thresh
                    t_old = t;
                    %w = A'*u;
                    t = K*u;
                    if norm(t)<thresh
                        flag = 1;
                        reason = 2;
                        break;
                    end

                    if itrns>15000
                        flag = 1;
                        reason = 4;
                        break;
                    end

                    s = norm(t);
                    t = t/s;

                    c = L'*t;
                    u_old = u;
                    u = L*c;
            %         if norm(u)<thresh
            %             flag = 1;
            %             break;
            %         end
                    u = u/norm(u);
                    itrns = itrns+1;
                end
                itrns
                if flag
                    dim = indx-1
                    break;
                end

                U = [U,u];
                T = [T,t];
                S = [S,s];
                C = [C,c];
                %w = A'*u_old;
                %w = w/norm(w);
                %W = [W,w];


                v = (eye(n)-t*t');
                K = v*K*v;
                %A = v*A;
                L = v*L;
                normK = norm(K);
                if normK<stp%thresh
            %         %flag = 1;
                    reason = 3;
                    dim = indx;
                    break;
                end
                indx = indx +1;
            end
            %% % get the last dimension - not good this diemnsion is not completely orthogonal to
            %%% others
            % u = rand(n,1);
            % t = K*u;
            % t = t/norm(t);
            % c = L'*t;
            % u = L*c;
            % u = u/norm(u);
            % t = K*u;
            % t = t/norm(t);
            % U = [U,u];
            % T = [T,t];
            % C = [C,c];
            %%

            disp 'Reason of exit loop is'
            reason

            disp 'Number of iteration:'
            dim

            MAP.T_ = T;
            MAP.U_ = U;
            
            if dim>0
                MAP.D_ = size(T,2);
                Y1 = MAP.K*MAP.M2;
            else
                MAP.D_ = 0;
                Y1 = zeros(size(X,1),1);
            end
            
        end
        %%%%%%% end learn function
        
        %% input output functions
        function obj = setDimension(obj,d)
            if d>0
                obj.D_ = d;
            else
                obj.D_ = size(obj.T_,2);
            end
        end
        
        function d = getDimension(obj)
             d = obj.D_;
        end
        
        function r = T(obj)
            if isempty(obj.T_)
                error('You have to learn first.');
            end
            
            dims = 1:obj.D_;
            r = obj.T_(:,dims);
        end
        
        function r=U(obj)
            if isempty(obj.U_)
                error('You have to learn first.');
            end
            dims = 1:obj.D_;
            r = obj.U_(:,dims);
        end
        
        function r=get.R2(obj)
            if isempty(obj.T_)
                error('You have to learn first.');
            end
            
            r = obj.U*inv(obj.T'*obj.K*obj.U);
        end
        function r=R(obj)
            if numel(size(obj.X))==2 %if A1 is not tensor
                r = obj.X'*obj.R2;
            end
        end
        function r=get.M2(obj)
            r = obj.R2*obj.T'*obj.Y;
        end
        function r=M(obj)
            if numel(size(obj.X))==2 %if A1 is not tensor
                r = obj.X'*obj.M2;
            end
        end
        
        function r=Orig_X(obj)
            r = bsxfun(@plus,obj.X, obj.mean_orig_X);
        end
        
        function r=Orig_Y(obj)
            r = bsxfun(@plus,obj.Y, obj.mean_orig_Y);
            %r = max(r,[],2);
        end
    end
end