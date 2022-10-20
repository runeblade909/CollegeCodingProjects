function [p,q,r] = rollrate2pqr(phi,theta,psi,phi_dot,theta_dot,psi_dot)

    %% Sovle for pqr vector
    pqr = [phi_dot,theta_dot,psi_dot]...
        *inv([1, sin(phi)*tan(theta), cos(phi)*tan(theta);...
        0, cos(phi), sin(phi);...
        0, sin(phi)*sec(theta), cos(phi)*sec(theta)]);

    %% Define p, q, and r based on previous calculations
    p = pqr(1);
    q = pqr(2);
    r = pqr(3);
    
end