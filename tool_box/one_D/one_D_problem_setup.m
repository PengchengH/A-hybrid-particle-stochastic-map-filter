%% ----------------------------
%% SETUP DA PROBLEM
%% ----------------------------
% Define model parameters
Cycle     = 100;
T_Steps   = 100;
T_BurnIn  = 1;
T_SpinUp  = 1;

d         = 1;

sigma_x   = 1;  % noise variance
Q         = sigma_x^2*eye(d);
sigma_y   = 2.5;
R         = sigma_y^2*eye(d);

dt        = 1;
ds_obs    = 1;
dt_iter   = 1;

lik       = @(x,y,sig) scalar_Gaussian(x,y,sig); % likelihood
loglik    = @(x,y,sig) log_scalar_Gaussian(x,y,sig); % log likelihood
sampleLik = @(x,sig) x(:) + sig*randn(length(x),1);  % sample likelihood

% set initial condition for data generation & spin-up
m0 = 20;
C0 = 1;

% Setup forward operator
for_op = @(xt, t, dt, K) one_D_signal(xt, t);

% Setup observation operator
obs_s  = 1:ds_obs:d;
nobs_s = length(obs_s);
H = eye(d); H = H(obs_s,:);
obs_op = @(x) H*x;

% Save model parameters
model = struct;
model.name      = 'one_D';
model.Cycle     = Cycle;
model.d         = d;
model.dt        = dt;
model.dt_iter   = dt_iter;
model.R         = R;
model.Q         = Q;
model.lik       = lik;
model.loglik    = loglik;
model.sampleLik = sampleLik;
model.data_idx  = obs_s;
model.m0        = m0;  % initial conditions
model.C0        = C0;

% save number of steps
model.T_BurnIn  = T_BurnIn;
model.T_Steps   = T_Steps;
model.T_SpinUp  = T_SpinUp;

% Save operators
model.for_op    = for_op;
model.ds_obs    = ds_obs;
model.obs_op    = obs_op;
model.obs_lin   = H;

% Save simulation parameters
model.seed      = sd;

fprintf('Setup DA problem\n')

%% ----------------------------
%% GENERATE DATA
%% ----------------------------
% set initial condition
x0 = zeros(d, Cycle);
for i = 1:Cycle
    x0(:,i) = (m0 + sqrtm(C0)*randn(d,1))';
end
model.x0 =x0;


% run dynamics and generate data
model = generate_data(model, T_SpinUp + T_Steps);

fprintf('Generated data\n')

% -- END OF FILE --