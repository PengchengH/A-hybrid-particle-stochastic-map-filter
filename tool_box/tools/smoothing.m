function [xf] = smoothing(xf,h)
   Q = cov(xf);
   K = size(xf,1);
   d = size(xf,2);
   [~,p]=chol(Q);     
   if p==0
      a   = sqrt(1-h^2);
      xf  = a*xf + (1-a)*mean(xf,1);
      xf  = xf   + h*mvnrnd(zeros(d,1),Q,K);     
   end
end

