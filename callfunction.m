function [] = callfunction(num_vry,model_store,method_store,...
                M_vect,theta_v,orderTM_vect,rho,mode)

%     num_vry    = [4 3 2 5];
    M_c        = num_vry(1);
    method_c   = num_vry(2);
    order_c    = num_vry(3);
    theta_c    = num_vry(4);
   
    
    option = struct;
    option.M       = M_vect(M_c);
    option.theta   = theta_v(theta_c);
    option.orderTM = orderTM_vect(order_c);
    option.rho     = rho;
    
    mena = method_store{method_c}.name;
    
%     if strcmp(mena,'EnKF') && M_c ==4

    if strcmp(mena,'PSMF') || strcmp(mena,'PSMF_nsm') ||...
       (strcmp(mena,'SMF')&& theta_c == 1) || ...
       ((strcmp(mena,'SIR_ESRF') || strcmp(mena,'GMM_EnKF')) && order_c == 1)||...
       (order_c == 1 && theta_c == 1)
   
            callfunction_further( model_store{M_c}, method_store{method_c},mode,option); 
            
    end
%     end

end

function [] = callfunction_further(model, method,varargin)

% set fixed parameters and the default values
fixpara = struct;
fixpara.distMat       = metric_lorenz(model.d); % Pairwise distance matrix
fixpara.nonId_radius  = model.d;  % localisaiton radius
fixpara.offdiag_rad   = model.d;  % map off diagial radius_parse map

% default value set
fixpara.M   = 500;
fixpara.rho = 0; 
if strcmp(method.name,'PSMF') || strcmp(method.name,'PSMF_nsm') || ...
   strcmp(method.name,'SIR_ESRF')|| strcmp(method.name,'GMM_EnKF') 
    fixpara.theta = 0.5;
end
if strcmp(method.name,'SMF') || strcmp(method.name,'PSMF')|| ...
   strcmp(method.name,'PSMF_nsm')  
    fixpara.order_all = 1;  % order of map
end

% extract paramters
if ~isempty(varargin)
    if isstruct(varargin{end})
        opti = varargin{end};
        
        if isfield(opti,'M')
            fixpara.M = opti.M;     
        end        
        if isfield(opti,'rho')
            fixpara.rho = opti.rho;         
        end
        
        if strcmp(method.name,'PSMF') || strcmp(method.name,'PSMF_nsm') ||...
           strcmp(method.name,'SIR_ESRF')|| strcmp(method.name,'GMM_EnKF') 
            if isfield(opti,'theta')
                fixpara.theta = opti.theta;
            end
        end       
        if strcmp(method.name,'SMF') || strcmp(method.name,'PSMF') ||...
           strcmp(method.name,'PSMF_nsm')    
            if isfield(opti,'orderTM')
                fixpara.order_all = opti.orderTM;  % order of map
            end
        end     
        
    else  
        error('option should be a struct');      
    end
end

% run filter for all parameters
if ~isempty(varargin)
    if ~isstruct(varargin{1})
        if strcmp(varargin{1},'concise') || strcmp(varargin{1},'full')
            filter = run_filter_params(model, method.filter,fixpara,varargin{1});
        else
            error('the mode should be concise or full');
        end
    else
        filter = run_filter_params(model, method.filter,fixpara);
    end
end

% save workspace
if strcmp(method.name,'PSMF') || strcmp(method.name,'PSMF_nsm')
    
    save(['result/' model.name '/' method.name '_M' num2str(fixpara.M)...
         '_theta' num2str(fixpara.theta) '_order' num2str(fixpara.order_all)...
         '.mat'], 'filter'); 
    fprintf([method.name ' M = %d, theta= %d, order= %d '],...
         fixpara.M, fixpara.theta, fixpara.order_all);
    
elseif strcmp(method.name,'SIR_ESRF')|| strcmp(method.name,'GMM_EnKF')
    
    save(['result/' model.name '/' method.name '_M' num2str(fixpara.M)...
         '_theta' num2str(fixpara.theta) ...
         '.mat'], 'filter');  
    fprintf([method.name ' M = %d, theta= %d'],...
         fixpara.M, fixpara.theta);   
     
elseif strcmp(method.name,'SMF')
    
    save(['result/' model.name '/' method.name '_M' num2str(fixpara.M)...
          '_order' num2str(fixpara.order_all)...
         '.mat'], 'filter'); 
    fprintf([method.name ' M = %d, order= %d '],...
         fixpara.M, fixpara.order_all);   
     
else
    
    save(['result/' model.name '/' method.name '_M' num2str(fixpara.M)...
         '.mat'], 'filter');   
    fprintf([method.name ' M = %d '], fixpara.M);  
    
end

end



