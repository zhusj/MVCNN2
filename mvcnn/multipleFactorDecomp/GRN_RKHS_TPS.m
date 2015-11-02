% learn Kernel PLS mapping: implementatin to what in paper Kenrk Partial
% Least square regression in RKHS.
% input:
% D : original data matrix Nxd  N: is the number of samples, d: is the dimensionality of the input images
% X : embedding space representation. Nxe matrix N: number of samples, e: the dimesnionality of the embedding space
% cent: centers: exm e: the dimesnionality of the embedding space

classdef GRN_RKHS_TPS
    % implementation of Generalized R... Network in paper []
    properties
        Cents=[];%Centers
        Map%map matrix -- regressor
        n_cents=0;
        Lambda
    end
    
    properties(SetAccess = private)
        %%%
    end
    properties(Dependent = true, SetAccess= private)
        %%%
    end
%     
    methods
%         function obj = GRN_RKHS(obj,p_cents,p_map)
%             obj.Cents = p_cents;
%             obj.Map = p_map;
%             obj.n_cents = size(p_cents,1);
%         end
        
        function obj = set_map(obj,p_map)
            obj.Map = p_map;
        end
        
        function obj = set_cents(obj,p_cents)
            obj.Cents = p_cents;
            obj.n_cents = size(p_cents,1);
        end
        
        
        function [Y]=get_kernel(obj,X)
            cents = obj.Cents;
            %X = p_input;
            Y =full(directs_grbf_reg(X,cents));
        end
        
        function [Y]=predict_regularized(obj,p_input)
            cents = obj.Cents;
            dim = size(cents,2);
            X = p_input;
            %implementation of higher dimensional regularization
            %cf = (G'G+lambdag)-1 G'y
            G=full(directs_grbf_reg_poly(X,cents,obj.Lambda));
            %G = bsxfun(@rdivide,G,sum(G,2));
            if numel(size(obj.Map))==2
                Y = G*obj.Map ;
            else
                Y = tmul(obj.Map,G',3);
            end
            Y = Y(1:end-dim-1,:);
            
        end
        
        function [obj,map]=learn_regularized(obj,p_input,p_output,p_n_cents,lambda,p_cents)
            % learn GRBF mapping
            % input:
            % output : original data matrix Nxd  N: is the number of samples, d: is the dimensionality of the input images
            % input :
            % n_cent: centers: exm e: the dimesnionality of the embedding space
            %I have input points, output points and need regression matrix.

            %% computing the centers of the RKHS
            if(exist('p_cents','var') && ~isempty(p_cents))
                obj.Cents = p_cents;
                obj.n_cents = size(obj.Cents,1);
            else
                if(~exist('p_n_cents','var') || p_n_cents<1)
                    obj.n_cents = 2*round(0.1*size(p_input,1));%this number of centers i
                    %if(obj.n_cents <8)
                    %    obj.n_cents = 8;
                    %end
                else
                    obj.n_cents = p_n_cents;
                end
                % centers are the centers of kmean of the input points.
                if obj.n_cents < size(p_input,1)
                    [~,C_] = kmeans(p_input,obj.n_cents);
                    obj.Cents = C_;
                else
                    obj.Cents = p_input;
                    obj.n_cents = size(obj.Cents,1);
                end
            end
%            %compute centers as the average of latent points for every class
%             unique_lbls = unique(p_labels);
%             cent = zeros(numel(unique_lbls),size(p_input,2));
%             c = 1;
%             for i = unique_lbls'
%                 indx = (p_labels==i);
%                 l = p_input(indx,:);
%                 cent(c,:) = mean(l);
%                 c = c+1;
%             end
            
            %% start the GRN regression
            D = p_output;
            cents = obj.Cents;
            X = p_input;
            %N=size(D,2);
            if(~exist('lambda','var'))
                lambda = -1.0;
            end
            obj.Lambda = lambda;
            %thresh = 1E-5;
            
%             %implementation of higher dimensional regularization
%             %cf = (G'G+lambdag)-1 G'y
%             G =directs_grbf_reg(X,cents);
%             %G = bsxfun(@rdivide,G,sum(G,2));
%             g =directs_grbf_reg(cents,cents);
%             %g = bsxfun(@rdivide,g,sum(g,2));
%             
%             F =G'*D;
%             B =full(G'*G-lambda*g);
%             map = B\D;
            [N,q] =size(p_input);
            d = size(D,2);
            G=full(directs_grbf_reg_poly(X,cents,lambda));
            
            D=[D; zeros(q+1,d)];
            
            obj.Map = G\D;
            
            %obj.n_cents
            
            %D = G*map;
            %disp 'The reconstruction error is:'
            %D = D(1:end-q-1,:);
            %err = norm(D - p_output);
            %obj.Map = map;
            %%%%%
            % CF = CF(2:end,:);
            %CF = inv(G'*G-lambda*g)*F;
        end
            
        
    end
end