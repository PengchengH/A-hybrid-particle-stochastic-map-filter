%% ----------------------------
%% SETUP DA PROBLEM
%% ----------------------------
% Define model parameters
Cycle     = 1;
T_Steps   = 1000;
T_BurnIn  = 1;
T_SpinUp  = 1;

x_real=[-4;.6];

d         = 2;

sigma_x = 0;  % noise variance
Q       = sigma_x^2*eye(d);
sigma_y = diag([1;0.1]);
R       = sigma_y^2*eye(d);

dt        = 1;
ds_obs    = 1;
dt_iter   = 1;

lik       = @(x,y,sig) scalar_Gaussian(x,y,sig); % likelihood
loglik    = @(x,y,sig) log_scalar_Gaussian(x,y,sig); % log likelihood
sampleLik = @(x,sig) x(:) + sig*randn(length(x),1);  % sample likelihood

% set initial condition for data generation & spin-up
m0 = zeros(d,1);
C0 = eye(d);

% Setup forward operator
for_op = @(xt, t, dt, K) Henon(randn(2,K))';

% Setup observation operator
obs_s  = 1:ds_obs:d;
nobs_s = length(obs_s);
H = eye(d); H = H(obs_s,:);
obs_op = @(x) H*x;

% Save model parameters
model = struct;
model.name      = 'Henon';
model.Cycle     = Cycle;
model.d         = d;
model.dt        = dt;
model.dt_iter   = dt_iter;
% model.sigma_x   = sigma_x;
% model.sigma_y   = sigma_y;
model.R         = R;
model.Q         = Q;
model.lik=lik;
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

model.x_real =x_real;
model.xt= diag(x_real)*ones(d,T_Steps+T_SpinUp);
model.yt=x_real+sigma_y*randn(d,T_Steps+T_SpinUp);

% Save simulation parameters
model.seed      = sd;

fprintf('Setup DA problem\n')

fprintf('Generated data\n')

% -- END OF FILE --