function filter_out = run_filter_params(model, method, fixed_params,varargin)

% create fixed parametes list
fn_list = fieldnames(fixed_params)';
nfn_list=size(fn_list,2);
    
% add fixed_params to options
for j=1:nfn_list
    options.(fn_list{j})=fixed_params.(fn_list{j});
end

% run filter and post-process outputs
algorithm = method(model, options);   % define class

if ~isempty(varargin)
    if strcmp(varargin{1},'concise')
        filter = seq_assimilation_concise(model, @algorithm.sample_posterior, options.M);
    elseif strcmp(varargin{1},'full')
        filter = seq_assimilation_full(model, @algorithm.sample_posterior, options.M);
    else
        error('the mode name should be concise or full');
    end
end

filter_out = post_process(model, filter, options); 

end
