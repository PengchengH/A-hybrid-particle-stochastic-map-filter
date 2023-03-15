function [y] = resampling(x,w)

K = size(x,1);
y = x;

F_w = zeros(1,K);
F_w(1) = w(1);
for i = 1:K-1
     F_w(i+1) = F_w(i)+ w(i+1);
end     

j = 1;          
for i=1:K  
    samp  = (i-1+rand)/K; 
    pp    = j;
    for j = pp:K
        if F_w(j)>=samp
            y(i,:)=x(j,:);
            break;
        end
    end
end

end

