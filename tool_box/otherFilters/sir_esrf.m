classdef  sir_esrf

	properties

		model      % struct containing definition of observations
		options    % struct containing options for EnKF algorithm

	end

	methods
		function obj =  sir_esrf(model, options, varargin)

			% define EnKF object
			p = ImprovedInputParser;
            addRequired(p,'model');
            addRequired(p,'options');
			parse(p, model, options, varargin{:});
			obj = passMatchedArgsToProperties(p, obj);

		end %endFunction
		% -----------------------------------------------------------------
		% -----------------------------------------------------------------
		function [out] = sample_posterior(obj, Xpr, Yt)

			% Load model inputs
			H     = obj.model.obs_lin;
            R     = obj.model.R;
            theta = obj.options.theta;
            K     = size(Xpr,1);

			% Load EnKF inputs
			rho       = obj.options.rho;

			% extract subset of data and data_idx
			data_idx  = obj.model.data_idx;
            n_obs     = length(data_idx);
            
            % looking for a proper alpha based on Ess
            alpha_start = 0;
            alpha_end   = 1;
            if obj.ESS_Calculation(Xpr, alpha_end, Yt) >= K*theta
               alpha = alpha_end;
            else
               for i = 1:10
                   alpha = 0.5*(alpha_start + alpha_end);
                   Ess   = obj.ESS_Calculation(Xpr, alpha, Yt);
                   if Ess >= K*theta
                      alpha_start = alpha;
                   else
                      alpha_end   = alpha;
                   end 
               end
               alpha = 0.5*(alpha_start + alpha_end);
            end
           
            % forcast of observation and calculate the weights
            Zpr = Xpr*H';
            w_p = zeros(K,1);     
            for i = 1:K       
                w_p(i) = exp( -0.5*(Zpr(i,:)-Yt')...
                         /(R/alpha)*(Zpr(i,:)-Yt')' ); %%weights
            end
            
            w_sum = sum(w_p);
            if w_sum == 0 || isnan(w_sum)
               w_p = ones(1,K)/K;
            else
               w_p = w_p/w_sum; 
            end
                 
            % resampling
            [Xpos] = resampling(Xpr,w_p);

			% run assimilation sequentially
			xf = Xpos;
            
            % inflaction
            mean_x = mean(xf,1);
            xf     = sqrt(1 + rho)*(xf - mean_x) + mean_x; 
            
            for i = 1:n_obs
                
				% extract local data and local observation operator
				data_local = Yt(i,:);
                 
                % Calculate mean;
                mean_x = mean(xf,1);
		        
				% sample from conditional likelihood using inflated particles               
                A   = ((xf - mean_x) / sqrt(K-1))';
                V   = A(i,:)';
                s2  = V'*V;   
                sigma_y = sqrtm(R(data_idx(i),data_idx(i)));
                sigma_y = sigma_y / sqrt(1-alpha);
                g2  = sigma_y^2;
                WHB = 1/(s2 + g2 + sqrt(s2*g2 + g2^2));
                
                A_a = A - WHB*A*V*(V');
                x_a = mean_x' + (data_local - mean_x(i))/(s2+g2)*A*V;
                xf = (x_a + sqrt(K-1)*A_a)';
                               
            end
                TH  = getMPRR(K);

			% set Xpost
            out  = struct;
            out.xpos  = TH'*xf;
            out.alpha = alpha;           

		end %endFunction
		% -----------------------------------------------------------------
        %------------------------------------------------------------------
        function [N_eff] = ESS_Calculation(obj, X, alpha, Yt)
            % Load model inputs
			H         = obj.model.obs_lin;
            R         = obj.model.R;
            K         = size(X,1); 
            
            Z = X*H';
            w = zeros(K,1);     
            for i = 1:K       
                w(i) = exp(-0.5*(Z(i,:)-Yt') / (R/alpha) * (Z(i,:)-Yt')'); %%weights
            end
            
            % Balance weights and resampling
            w_sum = sum(w);
            if w_sum == 0 || isnan(w_sum)
               w     = ones(1,K)/K;
               w_sum = 1;
            end
            
            % effective number of particles
            N_eff = 0;
            for i = 1:K
                w(i)  = w(i)/w_sum; 
                N_eff = N_eff+w(i)^2;
            end
            N_eff = 1/N_eff;       
%             plot(w);           
        end
        %------------------------------------------------------------------      
		% -----------------------------------------------------------------
	end %endMethods
end %endClass