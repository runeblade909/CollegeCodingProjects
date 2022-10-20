function [x,t] = linearizedEOM(x0,t_f)

    g = 9.81; %acceleration due to gravity
    
    %% Defien state space matrix
    A = zeros(12,12);
    A(1,7) = 1;
    A(2,8) = 1;
    A(3,9) = 1;
    A(4,10) = 1;
    A(5,11) = 1;
    A(6,12) = 1;
    A(7,5) = -g;
    A(8,4) = g;
    
    %% Define state space system
    sys = ss(A,zeros(12,1),eye(12),0);
    
    %% Simulate state space system
    [x,t] = initial(sys,x0,t_f);    

end