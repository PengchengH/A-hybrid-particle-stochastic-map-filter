clear; clc; close all; warning('off');

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

% state-space models
SS_model = {'L63', 'Henon', 'one_D','TT1'};

for SS_model_c = 1 : size(SS_model,2)

    % Build initial conditions      
    if strcmp(SS_model{1,SS_model_c},'L63')
        L63_problem_setup;
    elseif strcmp(SS_model{1,SS_model_c},'Henon')
        Henon_problem_setup;
    elseif strcmp(SS_model{1,SS_model_c},'one_D')
        one_D_problem_setup;
    elseif strcmp(SS_model{1,SS_model_c},'TT1')
        TT_problem_setup;  
    end

    model_init   = model;

    % define model save matrix;
    model_store = cell(size(M_vect,2),1);

    % Run spin-up and save data
    for M_c = 1:size(M_vect,2)   
        % produce prior samples by EnKF;
        M      = M_vect(M_c); 
        model  = spin_up(model_init, M);  

        path = ['result/' SS_model{1,SS_model_c}];
        mkdir (path);

        % save files
        save( [ 'result/' SS_model{1,SS_model_c}  ...
              '/spinup_M' num2str(M)], 'model');
        model_store{M_c,1} = model;
    end
        disp('Spin-up complete!')
end