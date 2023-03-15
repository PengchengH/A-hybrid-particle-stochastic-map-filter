clear; clc;  warning('off');

sd = 10; rng(sd);

% add functions to current path
addpath('tool_box/dataAssimilation')
addpath('tool_box/methods')
addpath('tool_box/otherFilters')
addpath('tool_box/tools')
addpath('tool_box/L63')
addpath('tool_box/Henon')
addpath('tool_box/one_D')
addpath('tool_box/TargetTracking')

% define number of samples
M_vect = [10,20,40,60,100,200,400,600];

% define tuning parameter for inflation factor
rho = 0;

% define alpha
theta_v =  [0.001,0.002,0.004,0.006,0.008, ...
           0.01, 0.02, 0.04, 0.06, 0.08,  ...
           0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,...
           0.92, 0.94, 0.96, 0.98, 0.99,  ...
           0.992,0.994,0.996,0.998,0.999];

% define orders for stochastic maps
orderTM_vect = [1,3];


% define number of processors
n_proc   = 8;

% method definition
method_store           = cell(8,1);
method_store{1}        = struct;
method_store{1}.filter = @pf;
method_store{1}.name   = 'PF';
method_store{2}        = struct;
method_store{2}.filter = @smf;
method_store{2}.name   = 'SMF';
method_store{3}        = struct;
method_store{3}.filter = @psmf;
method_store{3}.name   = 'PSMF';
method_store{4}        = struct;
method_store{4}.filter = @enkf;
method_store{4}.name   = 'EnKF';
method_store{5}        = struct;
method_store{5}.filter = @esrf;
method_store{5}.name   = 'ESRF';
method_store{6}        = struct;
method_store{6}.filter = @sir_esrf;
method_store{6}.name   = 'SIR_ESRF';
method_store{7}        = struct;
method_store{7}.filter = @gmm_enkf;
method_store{7}.name   = 'GMM_EnKF';
method_store{8}.filter = @psmf_nsm;
method_store{8}.name   = 'PSMF_nsm';


% state-space models
SS_model = {'L63', 'Henon', 'one_D','TT1'};

for SS_model_c = 1 : size(SS_model,2)

    % define model save matrix;
    model_store = cell(size(M_vect,2),1);

    % Load data
    for M_c = 1:size(M_vect,2)   
        M      = M_vect(M_c); 
        % save files
        load( [ 'result/' SS_model{1,SS_model_c} ...
              '/spinup_M' num2str(M)], 'model');
        model_store{M_c} = model;
    end
    
    disp('load data complete!');

    % variable number set
    num_base = [size(M_vect,2), size(method_store,1),...
               size(orderTM_vect,2), size(theta_v,2)];
    num_vry_index = parallel_arrange(num_base);

    % Parallel computation
    c = parcluster('local');   % set clusters
    c.NumWorkers = n_proc;
    parpool(n_proc);
    parfor (i_index = 1:size(num_vry_index,1), n_proc)
%     for i_index = 1:size(num_vry_index,1)

        callfunction(num_vry_index(i_index,:),model_store,method_store,...
             M_vect,theta_v,orderTM_vect,rho,'full');
 
    end

    % delete parallel pool object
    delete(gcp('nocreate'))

end
% -- END OF FILE --
