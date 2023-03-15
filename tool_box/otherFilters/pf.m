classdef pf
	% classic particle filter

	properties

		model      % struct containing definition of observations
		option    % struct containing option for EnKF algorithm

	end

	methods
        function obj = pf(model, option, varargin)
            
            % check the input parameter
			p = ImprovedInputParser;
            addRequired(p,'model');
            addRequired(p,'option');
			parse(p, model, option, varargin{:});
			obj = passMatchedArgsToProperties(p, obj);
            
            % default values of parameters
            % smoothing factor
            obj.option.smoothfactor = 0.2; 

		end %endFunction
		

		function [out] = sample_posterior(obj, xpr, yt)
            
			% Load model inputs
			obs_lin   = obj.model.obs_lin;
            R         = obj.model.R;
          
            % forcast of observation
            zpr = xpr*obs_lin';
            K   = size(xpr,1);
            wp = zeros(K,1);
            
            % weights
            for i=1:K       
                wp(i)=exp(-0.5*(zpr(i,:)-yt')/R*(zpr(i,:)-yt')'); 
            end
           
            % normalize weights   
            w_sum=sum(wp);
            if w_sum==0
                wp = ones(1,K)/K;
            else
                wp = wp/w_sum; 
            end
            
            % resampling
            xf  = resampling(xpr,wp); 
            
            % smoothing
            xpos = smoothing(xf,obj.option.smoothfactor);
            
			% set out
            out  = struct;
            out.xpos = xpos;

		end %endFunction
		% -----------------------------------------------------------------
		% -----------------------------------------------------------------
	end %endMethods
end %endClass



