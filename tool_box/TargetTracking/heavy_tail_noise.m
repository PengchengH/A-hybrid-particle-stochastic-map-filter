function u = heavy_tail_noise(Q1,Q2,p1,p2,K)

    d   = size(Q1,1);
    u   = zeros(K,d);
    
    for i=1:K
    s_d = rand(1);
    if s_d<p1
       u(i,:) = mvnrnd(zeros(d,1),Q1);
    else
       u(i,:) = mvnrnd(zeros(d,1),Q2);  
    end 
    end
end
