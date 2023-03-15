function filter_out = post_process(model, filter, options)
% declare struct to save results
filter_out = struct;

% extract indices for assimilation cycles
n0 = round(model.t0/(model.dt*model.dt_iter)) + 1;
Acycles = n0:n0+model.J-1;

% extract model information
xt = model.xt(:,Acycles,:);

% set T_BurnIn
if isfield(model,'T_BurnIn')
    T_BurnIn = model.T_BurnIn;
else
    T_BurnIn = 0;
end

% compute root mean square error
xt_mean = filter.mean(:,T_BurnIn+1:end,:); 
xt_true = xt(:,T_BurnIn+1:end,:);

if isfield(filter,'CRPS')
   CRPS  = filter.CRPS(:,T_BurnIn+1:end,:);
end

filter_out.timer = filter.timer;

% estimation RMSE
if strcmp(model.name,'Henon')
    
    filter_out.RMSE(1) = sqrt( mean( sum((xt_mean(1,:,:) - xt_true(1,:,:)).^2,1), 'all') );
    filter_out.RMSE(2) = sqrt( mean( sum((xt_mean(2,:,:) - xt_true(2,:,:)).^2,1), 'all') );
    
elseif strcmp(model.name,'TT1')
    
    filter_out.RMSE  = sqrt( mean( sum((xt_mean(1:2,:,:) - xt_true(1:2,:,:)).^2,1), 'all') );  
    filter_out.RMSE_v = sqrt( mean (sum((xt_mean(3,:,:) - xt_true(3,:,:)).^2,1), 'all') ); 
    filter_out.RMSE_algle = sqrt( mean( sum((xt_mean(4,:,:) - xt_true(4,:,:)).^2,1), 'all' ) );
  
else
    
    filter_out.RMSE  = sqrt( mean( sum((xt_mean - xt_true).^2,1), 'all' ) );
    
end


% estimation CRPS
if isfield(filter,'CRPS')
    if strcmp(model.name,'Henon')    
        filter_out.CRPS(1) = mean(CRPS(1,:,:),'all');
        filter_out.CRPS(2) = mean(CRPS(2,:,:),'all');
    elseif strcmp(model.name,'TT1')
        filter_out.CRPS = mean(CRPS(1:2,:,:),'all');
    else
        filter_out.CRPS = mean(CRPS,'all');
    end
end


% copy options to filter_out
for fn = fieldnames(options)'
    filter_out.(fn{1}) = options.(fn{1});
end

end

% -- END OF FILE --
