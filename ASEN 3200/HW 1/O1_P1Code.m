%% Part a)
% Initial Conditions
r = [7642,170,2186]'*1000; %[km]
r_dot = [0.32,6.91,4.29]'*(1000); %[km/s]

% Constants
M_E = 5.972*10^24;
G = 6.674*10^-11;
mu = G*M_E;
%Time Span
t0 = 0;
tf = 13000;% Single Period
IC = [r;r_dot];

% Ode45 Function Call:
[t,s] = ode45(@ode_function,[t0 tf],IC);
figure();
plot3(s(:,1),s(:,2),s(:,3));

function s_dot = ode_function(t,s)
% Assigning inputs
x = s(1);
y = s(2);
z = s(3);
vx =s(4);
vy =s(5);
vz = s(6);
r = [x;y;z];
% Constants
M_e = 5.972*10^24;
G = 6.674*10^-11;
k = (norm(r))^3;
% Calculating r_dot_dot
r_dot_dot = -G*M_e/k * r;
ax = r_dot_dot(1);
ay = r_dot_dot(2);
az = r_dot_dot(3);
% Final Outptu
s_dot = [s(4);s(5);s(6);ax;ay;az];
end