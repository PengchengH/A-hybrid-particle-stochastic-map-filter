% Author: Ricardo Baptista and Alessio Spantini and Youssef Marzouk
% Date:   May 2019
%
% See LICENSE.md for copyright information
%

classdef enkf
	% Define EnKF object to Assimilate observation at time t by 
	% generating analysis samples from the prior samples using 
	% the sequential perturbed observation EnKF. In this EnKF,
	% the inflation factor is applied to the prior samples before
	% sampling from the likelihood and the covariances are computed
	% using inflated prior and (not-inflated) likelihood samples.
	%
	% Notes: - Ensemble KF requires linear observation operator

	properties

		model      % struct containing definition of observations
		options    % struct containing options for EnKF algorithm

	end

	methods
		function obj = enkf(model, options, varargin)

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
			sampleLik = obj.model.sampleLik;
            R         = obj.model.R;

			% Load EnKF inputs
			rho       = obj.options.rho;

			% extract subset of data and data_idx
			data_idx  = obj.model.data_idx;

			% Determine the number of observations
			n_obs = size(obs_lin,1);

			% run assimilation sequentially
			xf = Xpr;
            
			for i=1:n_obs

				% extract local data and local observation operator
				obs_lin_local = obs_lin(i,:);
				data_local = Yt(i,:);

				% sample from conditional likelihood using inflated particles
				mean_x   = mean(xf,1);
				x_infl   = sqrt(1 + rho)*(xf - mean_x) + mean_x;
                sigma_y  = sqrtm(R(data_idx(i),data_idx(i)));
				y_sample = sampleLik(obs_lin_local*x_infl',sigma_y);
				
				% join x,y samples
				xy_samples = [y_sample, x_infl];
				mean_xy    = mean(xy_samples,1);
				delta_xy   = xy_samples - mean_xy;
				
				% estimate covariances
				CovAll = 1/(obj.options.M-1)*(delta_xy'*delta_xy);
				CovXY  = CovAll(2:end,1);
				VarY   = CovAll(1,1);

				% compute Kalman gain
				K = CovXY/VarY;
			    
				% sample from conditional likelihood using un-inflated particles
				y_sample = sampleLik(obs_lin_local*xf',sigma_y);

				% Find analysis samples
				xf = (xf' - K*(y_sample' - data_local))'; 

			end

			% set out
            out = struct;
			out.xpos = xf;

		end %endFunction
		% -----------------------------------------------------------------
		% -----------------------------------------------------------------
	end %endMethods
end %endClass