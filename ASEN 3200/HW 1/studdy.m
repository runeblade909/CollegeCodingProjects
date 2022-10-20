%Attitude

syms phi theta;


A = [1 0 0; 0 cos(phi) sin(phi);0 -sin(phi) cos(phi)];

B = [cos(theta) 0 -sin(theta); 0 1 0; sin(theta) 0 cos(theta)];

c = A*B;