function u = Transition_1(t,u, u1,T) 
    u(:,1) = u(:,1)+ T*cos(u(:,4)).*u(:,3);
    u(:,2) = u(:,2)+ T*sin(u(:,4)).*u(:,3);
    u      = u -u1;
end
