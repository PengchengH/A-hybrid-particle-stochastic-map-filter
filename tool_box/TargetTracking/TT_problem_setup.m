%% Define state-space model
Cycle     = 50;
T_Steps   = 120;
T_BurnIn  = 1;
T_SpinUp  = 1;


dt_iter   = 1;
ds_obs    = 1;

d         = 4;   % dimension of state
d_m       = 2;   % dimension of measurement


%% processing noise
dt  = 1;

%% select the tracking model
C0 = diag([10^2 10^2 3^2 (pi/10)^2]);
m0= [0 0 30 0]';

px_svar = 0.1;      % variance of position
py_svar = 0.1; 
v_svar  = 0.1;       % variance of velocity
D_svar  = 2*pi/180;  % variance of direction

sigma_x_1   = diag([px_svar,py_svar,v_svar,D_svar]);  % noise variance
sigma_x_2   = diag([10,10,10,30])*sigma_x_1;
Q_1         = sigma_x_1^2;
Q_2         = sigma_x_2^2;
u3 = [0 0 0 0];

p1=0.85;
p2=0.15;

% Setup forward operator
for_op = @(xt, t, dt, K) Transition_1(t, xt, u3, dt)+heavy_tail_noise(Q_1,Q_2,p1,p2,K);

Q=zeros(4,4);

%% measurement noise variance
sigma_y   = 3;
R         = diag([sigma_y^2,sigma_y^2]);

lik       = @(x,y,sig) scalar_Gaussian(x,y,sig); % likelihood
loglik    = @(x,y,sig) log_scalar_Gaussian(x,y,sig); % log likelihood
sampleLik = @(x,sig) x(:) + sig*randn(length(x),1);  % sample likelihood


% Setup observation operator
obs_s  = 1:ds_obs:d_m;
nobs_s = length(obs_s);
H = eye(d); H = H(obs_s,:);
obs_op = @(x) H*x;

% Save model parameters
model = struct;
model.name      = 'TT1';
model.Cycle     = Cycle;
model.d         = d;
model.d_m       = d_m;
model.dt        = dt;
model.dt_iter   = dt_iter;
model.Q         = Q;
model.R         = R;

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



% figure(1);
% subplot(2,2,1);
% plot(model.xt(1,:,1),model.xt(2,:,1),'-*');
% title('position(m)');
% subplot(2,2,2);
% plot(model.xt(3,:,1),'-*');
% title('velocityx(m/s)');
% subplot(2,2,3);
% plot(model.xt(4,:,1),'-*');
% title('direction(rad)');
 
fprintf('Generated data\n')

% -- END OF FILE --