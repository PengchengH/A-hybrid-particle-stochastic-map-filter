% Author: Ricardo Baptista and Alessio Spantini and Youssef Marzouk
% Date:   May 2019
%
% See LICENSE.md for copyright information
%

%% ----------------------------
%% SETUP DA PROBLEM
%% ----------------------------
% Define model parameters
Cycle     = 1;
T_Steps   = 4000;
T_BurnIn  = 2000;
T_SpinUp  = 2000;

d         = 3;
beta      = 8/3;
rho       = 28;
sigma     = 10;

sigma_x   = 0.01;  % noise variance
Q         = sigma_x^2*eye(d);
sigma_y   = 2;
R         = sigma_y^2*eye(d);

dt        = 0.05;
ds_obs    = 1;
dt_iter    = 2*5;

lik       = @(x,y,sig) scalar_Gaussian(x,y,sig); % likelihood
loglik    = @(x,y,sig) log_scalar_Gaussian(x,y,sig); % log likelihood
sampleLik = @(x,sig) x(:) + sig*randn(length(x),1);  % sample likelihood

% set initial condition for data generation & spin-up
m0 = zeros(d,1);
C0 = eye(d);

% Setup forward operator
for_op = @(xt, t, dt, K) rk4(@(t,u,dt) lorenz63(t,u,dt,sigma,rho,beta), xt, t, dt);

% Setup observation operator
obs_s  = 1:ds_obs:d;
nobs_s = length(obs_s);
H = eye(d); H = H(obs_s,:);
obs_op = @(x) H*x;

% Save model parameters
model = struct;
model.name      = 'L63';
model.Cycle     = Cycle;
model.d         = d;
model.dt        = dt;
model.dt_iter   = dt_iter;
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