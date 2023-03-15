function filter = seq_assimilation_concise(model, analysis_alg, M)
% SEQ_ASSIMILATION: Function assimilates observations sequentially
% using the algorithm specified with the function handle in analysis_alg

% Load model inputs
for_op  = model.for_op;
Q       = model.Q;
d       = model.d;
J       = model.J;
dt      = model.dt;
t0      = model.t0;
x0      = model.x0;
Cycle   = model.Cycle;
% xt      = model.xt; 

dt_iter = model.dt_iter;

% Declare vectors to store results
mean_t 	  = zeros(d,J,Cycle);
time_t 	  = zeros(J,Cycle);

xp  = zeros(M,d,Cycle);

tic; % timer

for Cycle_count = 1:Cycle
    
    % Initialisation
    xf = x0(:,:,Cycle_count);
    tf = t0;
    n0 = round(t0/(model.dt*model.dt_iter)) + 1;
    Acycles = n0:n0+J-1;

    % Assimilation
    for n=1:length(Acycles)
    
        for i = 1:dt_iter
            xf = for_op(xf, tf + dt*(i-1), dt, M);
        end   
        
	    xf_rand = xf + randn(size(xf))*sqrtm(Q);  % add process noise
        
	    tf   = tf + dt*dt_iter;   
    
        % extract data at assimilation step n
        Yt   = model.yt(:,Acycles(n),Cycle_count);

	    % assimisation
	    out  = analysis_alg(xf_rand, Yt);  % this is the function called
        
        xp(:,:,Cycle_count) = out.xpos;
        
        
        % set next forecast samples
	    xf  = xp(:,:,Cycle_count);
        
        mean_t(:,n,Cycle_count) = mean(xp(:,:,Cycle_count),1);
	    time_t(n,Cycle_count)   = tf;
        
        
%         if length(Acycles) == model.T_Steps
%             if n>50 && Cycle_count>0
%                 if size(xf,2) == 1
%                     disp_1D (xf_rand, Yt, xp(:,:,Cycle_count), model, mean_t, n, Cycle_count, Acycles);
%                 elseif size(xf,2) == 2
%                     disp_2D (xf_rand, Yt, xp(:,:,Cycle_count), model, mean_t, n, Cycle_count, Acycles); 
%                 elseif size(xf,2) == 3
%                     disp_3D (xf_rand, Yt, xp(:,:,Cycle_count), model, mean_t, n, Cycle_count, Acycles);
%                 elseif size(xf,2) == 4
%                     disp_4D (xf_rand, Yt, xp(:,:,Cycle_count), model, mean_t, n, Cycle_count, Acycles); 
%                 end  
% 
%            end
%         end


    end
end

T_comsume = toc;
% log(T_comsume)

% Return results
filter = struct;
filter.timer     = T_comsume;
filter.xp        = xp;
filter.time      = time_t;
filter.mean      = mean_t;

end

% -- END OF FILE --