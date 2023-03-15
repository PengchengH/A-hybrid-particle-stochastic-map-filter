% Author: Ricardo Baptista and Alessio Spantini and Youssef Marzouk
% Date:   May 2019
%
% See LICENSE.md for copyright information
%

function model = generate_data(model, J)

% Load model inputs
d 		  = model.d;
for_op    = model.for_op;
Q         = model.Q;
R         = model.R;
dt        = model.dt;
dt_iter   = model.dt_iter;
sampleLik = model.sampleLik;
Cycle     = model.Cycle;
x0        = model.x0;

% find number of observations in space from linear operator
nobs = length(model.data_idx);

% Declare vectors to store state, observations and time
xt = zeros(d, J, Cycle);
yt = zeros(nobs, J, Cycle);
tt = zeros(J, 1, Cycle);

for Cycle_count = 1:Cycle
    % Initialize xf and tf
    xf = x0(:,Cycle_count)';
    tf = 0;

    % Generate true data
    for n=1:J

	    % run dynamics and save results
        for i=1:dt_iter
            xf = for_op(xf, tf + dt*(i-1), dt, 1);
        end
	    xf = xf + randn(size(xf))*sqrtm(Q);
	    tf = tf + dt*dt_iter;

	    % collect observations
	    yt(:, n, Cycle_count) = sampleLik(model.obs_lin*xf', sqrtm(R));

	    % Save results in vector
	    xt(:, n, Cycle_count) = xf;
	    tt(n, Cycle_count)    = tf;
        
    end
    
end

% Save data in model
model.xt = xt;
model.yt = yt;
model.tt = tt;

% plot(xt(:,:,1));

end

% -- END OF FILE --