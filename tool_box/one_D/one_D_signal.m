function y = one_D_signal(x, t)
     y = 0.5*x + 25*x./(1+x.^2) + 8*cos(1.2*t);
end
