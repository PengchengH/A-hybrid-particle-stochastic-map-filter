classdef gmm_enkf
	% Define GMM_EnKF object to Assimilate observation at time t by 
	% generating analysis samples from the prior samples using 
	% the sequential perturbed observation GMM_EnKF. In this GMM_EnKF,
	% the inflation factor is applied to the prior samples before
	% sampling from the likelihood and the covariances are computed
	% using inflated prior and (not-inflated) likelihood samples.
	%
	% Notes: - Ensemble KF requires linear observation operator

	properties

		model      % struct containing definition of observations
		options    % struct containing options for GMM_EnKF algorithm

	end

	methods
		function obj = gmm_enkf(model, options, varargin)

			% define GMM_EnKF object
			p = ImprovedInputParser;
            addRequired(p,'model');
            addRequired(p,'options');
			parse(p, model, options, varargin{:});
			obj = passMatchedArgsToProperties(p, obj);

		end %endFunction
		% -----------------------------------------------------------------
		% -----------------------------------------------------------------
		function [out] = sample_posterior(obj, xpr, yt)

			% Load model inputs
			H         = obj.model.obs_lin;
            d         = obj.model.d;
            n_obs     = size(H,1);
            R         = obj.model.R;
            theta     = obj.options.theta;
            M         = size(xpr,1);

			% Load GMM_EnKF inputs
			rho       = obj.options.rho;
            
            % looking for a proper alpha based on Ess
            alpha_start = 0;
            alpha_end   = 1;
            if obj.ESS_Calculation_2(xpr, alpha_start, yt) >= M*theta
               alpha = alpha_start;
            else
               for i = 1:10
                   alpha = 0.5*(alpha_start + alpha_end);
                   Ess   = obj.ESS_Calculation_1(xpr, alpha, yt);
                   if Ess >= M*theta
                      alpha_end = alpha;
                   else
                      alpha_start   = alpha;
                   end 
               end
               alpha = 0.5*(alpha_start + alpha_end);
            end
            
            if alpha>0
               % inflaction
               mean_x = mean(xpr,1);
			   x_inf  = sqrt(1 + rho)*(xpr - mean_x) + mean_x; 
            
               % gmm-enkf assimilation;
               R_alpha  = R/alpha;
               Cov_x    = cov(x_inf);
               K_alpha  = Cov_x*H'/(H*Cov_x*H'+R_alpha);
               Cov_enkf = (1/alpha)*K_alpha*R*K_alpha';
               K_end      = Cov_enkf*H'/(H*Cov_enkf*H' + R/(1-alpha));
             
               % mean of Gaussain Mixture and expanded likelihood
               Cov_alpha = H*Cov_enkf*H'+R/(1-alpha);
               nu_alpha  = xpr+ (yt'-xpr*H')*K_alpha';
               nu_end    = nu_alpha + (yt'-nu_alpha*H')*K_end';
            
               % Calculate weights
               w = zeros(M,1);     
               for i = 1:M       
                   w(i) = exp(-0.5*(nu_alpha(i,:)*H'-yt') / Cov_alpha...
                          * (nu_alpha(i,:)*H'-yt')' ); %%weights
               end
 
               w_sum = sum(w);
               if w_sum==0 || isnan(w_sum)
                   w = ones(1,M)/M;
               else
                   w = w/w_sum; 
               end
            
               % resampling
               F_w = zeros(1,M);
               F_w(1) = w(1);
               for i = 1:M-1
                   F_w(i+1) = F_w(i)+ w(i+1);
               end 
               
               xpost = zeros(M,d);
               j = 1;
               for i = 1:M  
                   samp  = (i-1+rand)/M;
                   pp    = j;
                   for j = pp:M        
                       if F_w(j) >= samp
                          xpost(i,:) = nu_end(j,:)+mvnrnd(zeros(n_obs,1),R_alpha,1)...
                              *K_alpha'*(eye(d)-K_end*H)'+...
                              mvnrnd(zeros(n_obs,1),R/(1-alpha),1)*K_end';
                          break;
                       end
                   end
               end
               
            
            else
               % forcast of observation
               zpr = xpr*H';
               K   = size(xpr,1);
               w_p = zeros(K,1);
            
               for i=1:K       
                  w_p(i)=exp(-0.5*(zpr(i,:)-yt')/R*(zpr(i,:)-yt')'); %%weights
               end
           
               %%%%%%%%%%%%%%%%%normalize weights   
               w_sum=sum(w_p);
               if w_sum==0
                   w_p=ones(1,K)/K;
               else
                   w_p =w_p/w_sum; 
               end
            
               % resampling
               [xf] = resampling(xpr,w_p);
               
               % smoothing;
               xpost = smoothing(xf,0.2);

            end
            
            % set out
            out  = struct;
            out.xpos  = xpost;
            out.alpha = alpha; 

		end %endFunction
		% -----------------------------------------------------------------      
        %------------------------------------------------------------------
        function [N_eff] = ESS_Calculation_1(obj, x, alpha, yt)
            % Load model inputs
			H   = obj.model.obs_lin;
            R   = obj.model.R;
            rho = obj.options.rho;
            M   = size(x,1); 
            
            % inflaction
            mean_x = mean(x,1);
			x_inf  = sqrt(1 + rho)*(x - mean_x) + mean_x; 
            
            % EnKF related matrix
            Cov_x    = cov(x_inf);
            K_alpha  = Cov_x*H'/(H*Cov_x*H'+R/alpha);
            Cov_enkf = (1/alpha)*K_alpha*R*K_alpha';
            
            % mean of Gaussain Mixture and expanded likelihood
            Cov_alpha = H*Cov_enkf*H'+R/(1-alpha);
            nu_alpha  = x+ (yt'-x*H')*K_alpha';
            
            % Calculate weights
            w = zeros(M,1);     
            for i = 1:M       
                w(i) = exp(-0.5*(nu_alpha(i,:)*H'-yt') / Cov_alpha...
                       * (nu_alpha(i,:)*H'-yt')' ); %%weights
            end
            
            %---- Balance weights and resampling
            w_sum = sum(w);
            if w_sum == 0 || isnan(w_sum)
               w     = ones(1,M)/M;
               w_sum = 1;
            end
            
            % effective number of particles
            N_eff = 0;
            for i = 1:M
                w(i)  = w(i)/w_sum; 
                N_eff = N_eff+w(i)^2;
            end
            N_eff = 1/N_eff;       
%           plot(w);           
        end
        %------------------------------------------------------------------  
        %------------------------------------------------------------------
        function [N_eff] = ESS_Calculation_2(obj, X, alpha, Yt)
            % Load model inputs
			H         = obj.model.obs_lin;
            R         = obj.model.R;
%           d         = SIR_ESRF.model.d;
            K         = size(X,1); 
            alpha     = 1-alpha;
            
            Z = X*H';
            w = zeros(K,1);     
            for i = 1:K       
                w(i) = exp(-0.5*(Z(i,:)-Yt') / (R/alpha) * (Z(i,:)-Yt')'); %%weights
            end
            
            %---- Balance weights and resampling
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