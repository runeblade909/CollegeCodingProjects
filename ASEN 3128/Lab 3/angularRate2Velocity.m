function [p,q,r] = angularRate2Velocity(phi,theta,psi,phi_dot,theta_dot,psi_dot)

    %% EOM 2 for angular motion
    A = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);...
        0 cos(phi) -sin(phi);...
        0 sin(phi)*sec(theta) cos(phi)*sec(theta)];
    out = A*[phi_dot;theta_dot;psi_dot];
    
    p = out(1);
    q = out(2);
    r = out(3);
    
end