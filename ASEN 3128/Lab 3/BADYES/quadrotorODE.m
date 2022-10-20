function xdot = quadrotorODE(t,state_vec,m,I_x,I_y,I_z,n,mu,Z_c,L_c,M_c,N_c)

    phi = state_vec(4); % roll angle
    theta = state_vec(5); % pitch angle
    psi = state_vec(6); % yaw angle
    u = state_vec(7); % x velocity
    v = state_vec(8); % y velocity
    w = state_vec(9); % z velocity
    p = state_vec(10); % roll velocity
    q = state_vec(11); % pitch velocity
    r = state_vec(12); % yaw velocity

    g = 9.81; %Gravity m/s^2
    M = -mu*norm([p;q;r])*[p;q;r]; % aerodynamic moment 
    V = norm([u;v;w]); % velocity
    F = -n*V*[u;v;w]; % aerodynamic force 

    MomCntl = [L_c; M_c; N_c]; % control moment 

    % transformation matrix 
    R = [cos(theta)*cos(psi) sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);...
        cos(theta)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);...
        -sin(theta) sin(phi)*cos(theta) cos(phi)*cos(theta)];

    % equations of motion 
    v_inertial = R * [u;v;w];
    euler_dot = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);...
        0 cos(phi) -sin(phi);...
        0 sin(phi)*(1/cos(theta)) cos(phi)*(1/cos(theta))]...
        * [p;q;r];
    a_inertial = [r*v-q*w; p*w-r*u; q*u-p*v]...
        + g*[-sin(theta);cos(theta)*sin(phi);cos(theta)*cos(phi)]...
        + (1/m)*F + (1/m)*[0;0;Z_c];
    a_angular = [((I_y-I_z)/I_x)*q*r; ((I_z-I_x)/I_y)*q*r; ((I_x-I_y)/I_z)*q*r]...
        + [(1/I_x)*M(1); (1/I_y)*M(2); (1/I_z)*M(3)]...
        + [(1/I_x)*MomCntl(1); (1/I_y)*MomCntl(2); (1/I_z)*MomCntl(3)];

    xdot = [v_inertial; euler_dot; a_inertial; a_angular]; % vecotr output for the integrals

end