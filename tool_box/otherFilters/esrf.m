classdef  esrf

	properties

		model      
		options 

	end

	methods
		function obj =  esrf(model, options, varargin)

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
			obs_lin   = obj.model.obs_lin;
            R         = obj.model.R;

			% Load EnKF inputs
			rho       = obj.options.rho;

			% extract subset of data and data_idx
			data_idx  = obj.model.data_idx;

			% Determine the number of observations
			n_obs = size(obs_lin,1);

			% run assimilation sequentially
			xf = Xpr;
            K  = size(Xpr,1);
            
            % inflaction
            mean_x = mean(xf,1);
            xf = sqrt(1 + rho)*(xf - mean_x) + mean_x; 
                            
            for i=1:n_obs

				% extract local data and local observation operator
				data_local = Yt(i,:);
                
                mean_x = mean(xf,1);
             
                A  = ((xf - mean_x) / sqrt(K-1))';
                V  = A(i,:)';
                s2 = V'*V;  
                
                sigma_y  = sqrtm(R(data_idx(i),data_idx(i)));
                g2       = sigma_y^2;
                WHB      = 1/(s2 + g2 + sqrt(s2*g2 + g2^2));
                
                A_a = A - WHB*A*V*(V');
                x_a = mean_x' + (data_local - mean_x(i))/(s2+g2)*A*V;
                xf  = (x_a + sqrt(K-1)*A_a)';
              
            end

			% set Xpost
            out  = struct;
            out.xpos = xf;

		end %endFunction
		% -----------------------------------------------------------------
		% -----------------------------------------------------------------
	end %endMethods
end %endClass