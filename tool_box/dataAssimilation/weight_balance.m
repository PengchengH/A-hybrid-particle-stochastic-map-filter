function w_out=weight_balance(w_in)
          
          K=size(w_in,1);
            %---- Balance weights and resampling
            w_sum = sum(w_in);
            if w_sum == 0
                 w_in = ones(1,K)/K;
                 w_sum = 1;
            end
            
            % effective number of particles
            N_eff = 0;
            for i = 1:K
                 w_in(i) = w_in(i)/w_sum; 
                 N_eff  = N_eff+w_in(i)^2;
            end
            N_eff=1/N_eff;
            
            % weights balance
            alpha=N_eff/K;     
            for i = 1:K
                 w_in(i) =alpha*w_in(i)+(1-alpha)/K; 
            end
            
            w_out=w_in;
end