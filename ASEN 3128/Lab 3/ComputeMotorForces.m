% ASEN 3128 Lab 3 
%
% function to calculate motor forces using control moments
%
%% 

function motor_forces = ComputeMotorForces(Z_c,L_c,M_c,N_c, R, km)

%motor_forces = [f1,f2,f3,f4]
motor_forces = [-1,-1,-1,-1;...
                -R/sqrt(2),-R/sqrt(2),R/sqrt(2),R/sqrt(2);...
                R/sqrt(2),-R/sqrt(2),-R/sqrt(2),R/sqrt(2);...
                km,-km,km,-km]\[Z_c;L_c;M_c;N_c];
            
end