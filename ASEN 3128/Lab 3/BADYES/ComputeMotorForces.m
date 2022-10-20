function motor_forces = ComputeMotorForces(Z_c,L_c,M_c,N_c,r,k_m)

    %% Moment to force matrix
    m2f = [-1, -1, -1, -1;...
        -r/sqrt(2), -r/sqrt(2), r/sqrt(2), r/sqrt(2);...
        r/sqrt(2), -r/sqrt(2), -r/sqrt(2), r/sqrt(2);...
        k_m, -k_m, k_m, -k_m];
    
    %% Vecotr of control moments and force
    f_vec = [Z_c, L_c, M_c, N_c];
    
    %% Solve for motor forces
    motor_forces = m2f\f_vec';

end