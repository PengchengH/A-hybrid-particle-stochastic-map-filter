% Author: Ricardo Baptista and Alessio Spantini and Youssef Marzouk
% Date:   May 2019
%
% See LICENSE.md for copyright information
%

classdef smf
	% Define an object from the composition of TransportMaps to
    % generate posterior (from the conditional distribution Y|X)
    % samples given samples from the prior and the likelihood model

	properties

		order      % order of transport maps
        model      % struct with the likelihood and observation data
        options    % struct with options for assimilation
		TM  	   % map from samples for approximating posterior
        NonIdComp  % list of non-identity components in map

    end

	methods 
		function obj = smf(model, options, varargin)

			% define obj object
			p = ImprovedInputParser;
            addRequired(p,'model');
            addRequired(p,'options');
			parse(p, model, options, varargin{:});
			obj = passMatchedArgsToProperties(p, obj);

            % define default parameters for obj
            obj = obj.default_tmap();

            % define component orders and identify non-identity components
            [obj.order, obj.NonIdComp] = obj.setup_tmap();
            
            % define TM object
            d = length(obj.order);
           
            obj.TM = TransportMap(d, obj.order, obj.options);  %%transport map design

            % set remaining components to be identity function
            IdComp = setdiff(1:obj.TM.d, obj.NonIdComp);
            for Ck=IdComp
                obj.TM.S{Ck} = obj.TM.S{Ck}.set_Id_comp();
            end
            
		end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function obj = default_tmap(obj)
        % function adds recommended options for stochastic map filter

            obj.options.scalingWidths    = 2.0; % scale
            obj.options.lambda           = 0.1; % initial 0
            obj.options.delta            = 0.1; % initial 1e-8
            obj.options.npoints_interp   = max(2000,2*obj.options.M); % number of points for interpolation
            obj.options.kappa            = 4;

            obj.options.data_order       = obj.options.order_all;
            obj.options.offdiag_order    = obj.options.order_all; % off diagnal map term order
            obj.options.diag_order_obs   = obj.options.order_all; % diagnal observe
            obj.options.diag_order_unobs = 1;  % diagnal unobserve
            obj.options.locLik           = 1;  % local likelihood

        end %endFunction
        %------------------------------------------------------------------
        %------------------------------------------------------------------
        function [order, NonIdComp] = setup_tmap(obj)
        % define the orders for the S map components based on the specified
        % options and model parameters {order and non-identity component index}
        % the map about observation can be directly set as identity map;

            % extract parameters from model and options
            d = obj.model.d; % state dimension
            distMat = obj.options.distMat; % pairwise distance matrix

            % set off-diagonal orders for the map
            offdiag_order = obj.options.offdiag_order;  % off diagnal order
            offdiag_rad   = obj.options.offdiag_rad; % map off diagial radius_parse map
            dist_to_order = [offdiag_order*ones(1,offdiag_rad), zeros(1,d-offdiag_rad)];

            % define cell to store orders
            order = cell(d+1,1);

            % extract ordering of variables based on distance
            [~, permutation] = sort(distMat(:,1));
            
            for i=1:d
                
                % for each variable, extract distance to observed node
                node_i = permutation(i);
                dist_i_1 = distMat(1,node_i);

                % define vector to store orders
                orders_i = zeros(1,i);

                % if dist_i_1 <= nonIdcomponents, component is identity
                if dist_i_1 <= obj.options.nonId_radius

                    % define orders for off-diagonal TM component
                    for j = 1:i-1

                        % for each previous variable, extract distance to node_i
                        node_j = permutation(j);
                        dist_i_j = distMat(node_i,node_j);

                        % compute order for variable based on distance
                        orders_i(j) = dist_to_order(dist_i_j);
                        
                    end

                    % define orders for diagonal TM component
                    if i == 1 % observed node
                        orders_i(i) = obj.options.diag_order_obs;
                    else % unobserved node
                        orders_i(i) = obj.options.diag_order_unobs;
                    end

                    % append order for dependence on data and save orders
                    if (i == 1 && obj.options.locLik == 1) || (obj.options.locLik == 0)
                        order{i+1} = [obj.options.data_order, orders_i];
                    else
                        order{i+1} = [0, orders_i];
                    end

                end
            end

            % define non-identity components of map
            NonIdComp = find(cellfun(@sum,order)~=0)';

            % set remaining components to order=[0,\dots,0,1]
            IdComp = setdiff(1:d+1, NonIdComp);
            for i=IdComp
                order{i} = [zeros(1,i-1), 1];
            end
            
        end %endFunction
		%------------------------------------------------------------------  
        %------------------------------------------------------------------
        function [out] = sample_posterior(obj, xpr, yt)

            % extract number of indices for assimilation
			data_idx = obj.model.data_idx;
            n_obs = length(data_idx);
                        
			% generate analysis ensemble by sequentially assimilating data
            for i = 1:n_obs
                xpr = obj.assimilate_scalar_obs(xpr, data_idx(i), yt(i));
            end
            
			% set out
            out  = struct;
            out.xpos = xpr;

        end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
		function Xpost = assimilate_scalar_obs(obj, Xpr, data_idx, Yt)
        % assimilate scalar observation of X(:,data_idx) based on the
        % observation in Yt and by permutating the samples

            % check if Yt is a scalar
            if ~isscalar(Yt)
                error('Observation should be a scalar in obj.evaluate')
            end
            
            R       = obj.model.R;
            sigma_y = sqrtm(R(data_idx,data_idx));

            %----- Permute the ensemble -----
            [~, permutation] = sort(obj.options.distMat(:,data_idx));
            Xpr = Xpr(:,permutation);
                        
            % apply multiplicative inflation to forecast samples
            Xinfl = obj.inflate(Xpr);

            % generate samples from local likelihood
            Xi = Xinfl(:,1);
            Yi = obj.model.sampleLik(Xi,sigma_y);

            % compute lower map components using separability
            obj.TM = obj.TM.optimize([Yi, Xinfl], obj.NonIdComp);
            
            % generate local-likelihood samples with un-inflated X samples
            Xi  = Xpr(:,1);
            Yi  = obj.model.sampleLik(Xi,sigma_y);
            
            % evaluate composed map T(y,x) at [Yi, Xpr] samples
            Xpost = obj.evaluate([Yi, Xpr], Yt);
                        
            %----- Inverse permutation -----
            Xpost(:,permutation) = Xpost;
            
%             Q=cov(Xpost);
%             [~,p]=chol(Q);     
%             if p==0
%                  Xpost=Xpost+0.1*mvnrnd(zeros(size(Xpost,2),1),Q,size(Xpost,1));
%             end
            
		end %endFunction
		%------------------------------------------------------------------
		
		%------------------------------------------------------------------
        function infl_X = inflate(obj, X)
        	mean_X  = mean(X,1);
            delta_X = sqrt(1 + obj.options.rho)*(X - mean_X);
            infl_X  = delta_X + mean_X;
        end %endFunction
		%------------------------------------------------------------------
		%------------------------------------------------------------------
        function Tx = evaluate(obj, YX, Yt)
        % evaluate the composed map  T(y,x) at the samples YX = [Y, X]

            % check if Yt is a scalar
            if ~isscalar(Yt)
                error('Observation should be a scalar in obj.evaluate')
            end

            % find total number of samples
            N = size(YX,1);
            d=obj.model.d;

            % evaluate forward map at YX samples
            EtaYX = obj.TM.eval_map(YX,[1,d+1]); % from prior to Gaussian

            % set observation and evaluate inverse map
            EtaYX(:,1) = repmat(Yt, N, 1); 
            YtX = obj.TM.eval_inv_map(EtaYX,[1,d+1]); % from Gaussian to posterior

            % extract posterior samples
            Tx = YtX(:,2:end);

        end %endFunction
        %------------------------------------------------------------------
		%------------------------------------------------------------------
		function dTx = gradient(obj, YX, Yt)
        % compute the gradient of the composed map T at the samples X using
        % chain rule: T(y,x) = S(y*,\cdot)^{-1} \circ S(y,x)

            % find total number of samples
            N = size(YX,1);
        
            % evaluate \nabla S(y,x) for YX inputs
            dSx = obj.TM.grad_x(YX);
            
            % evaluate \nabla S(y^{*},x) for (Y^{*},T(Y,X)) inputs
            Tyx = obj.evaluate(YX,Yt);
            dSx_ysT = obj.TM.grad_x([Yt*ones(N,1), Tyx]);
            
            % evaluate \nabla_{x} S^{-1}(x)|S using inverse function 
            % theorem as inverse of (\nabla_{x} S(y^{*},T(y,x)))
            inv_dSinvx = dSx_ysT(:,2:end,2:end);
            
            % evaluate gradient of T by composition
            dTx = zeros(N, obj.TM.d-1, obj.TM.d);
            for i=1:N
                dTx(i,:,:) = squeeze(inv_dSinvx(i,:,:))\squeeze(dSx(i,2:end,:));
            end
            
        end %endFunction
        %------------------------------------------------------------------
		%------------------------------------------------------------------
	end %endMethods
end %endClass