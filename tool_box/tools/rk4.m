function [u_n] = rk4(func, u, t, dt)
        % Compute intermediary values for k
        k1 = func(t, u, dt);
        k2 = func(t + dt/2, u + dt/2*k1, dt);
        k3 = func(t + dt/2, u + dt/2*k2, dt);
        k4 = func(t + dt, u + dt*k3, dt);
        % Compute updated values for u and t
        u_n = u + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end
