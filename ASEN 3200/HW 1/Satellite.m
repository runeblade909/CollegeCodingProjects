function dstatedt = Satellite(t,state)
%   S0 = [r r_dot];

% x = state(1);
% y = state(2);
% z = state(3);
% xdot = state(4);
% ydot = state(5);
% zdot = state(6);

M = 5.97219e24;
G = 6.67408e-11;
mu = G*M;


r = state(1:3);
rnorm = norm(r);
v = state(4:6);

vel = v;
accel = -mu/(rnorm^3).*r;


dstatedt = [vel;accel];



end