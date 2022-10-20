%% ASEN 3128 Lab 2 Description: Quadrotor EOM
% Code Authors: Zak Reichenbach,Ketan Kamat,Jack Iribarren,Brian Byrne
% Date of Submission: 9/21/2021
% Code Description: Simulating motion of a quadrotor under different
% conditions
%% Cleaning up after each run
clc; 
clear; 
close all; 
%% Constants
% Declaring global variables for access within functions
global m radius km Ix Iy Iz vCoef mu g f1 f2 f3 f4 Lc Mc Nc

g       = 9.81;         %[m/s^2]
m       = 0.068;        %[kg]
radius  = 0.06;         %[m]
km      = 0.0024;       %[N*m/N]
Ix      = 6.8*10^(-5);  %[kg*m^2]
Iy      = 9.2*10^(-5);  %[kg*m^2]
Iz      = 1.35*10^(-4); %[kg*m^2]
vCoef   = 1*10^(-3);    %[N/(m/s)^2]
mu      = 2*10^(-6);    %[N*m/(rad/s)^2]

%% Problem (1)

% Initial Conditions
pos_0 = [0,0,-5];                   % Position [m]
ori_0 = [0,0,0];                    % Orientation [rad]
vel_0 = [0,0,0];                    % Velocity [m/s]
omega_0 = [0,0,0];                  % Rates [rad/s]
s_0_1 = [pos_0 ori_0 vel_0 omega_0];

% Force/Moment Calculations
F = ((m*g)/4); %[N]
f1 = F; %[N]
f2 = F; %[N]
f3 = F;%[N]
f4 = F;%[N]
[Lc,Mc,Nc] = controlLMNc(radius,f1,f2,f3,f4,km); % Calculating  control moments: %[N m]

% ODE 45 Function Call
t_span = [0 10]; % Time Span [s]
[t1,s1] = ode45(@odefun_level,t_span,s_0_1);  

% Plots
plot12(t1,s1);
fprintf("Problem 1: The trim thrust values are %d [N]\n",F);
%% Problem (2a)

% Initial Conditions
pos_0 = [0,0,-5];                   % Position
ori_0 = [0,0,0];                    % Orientation
vel_0 = [0,0,0];                    % Velocity
omega_0 = [0,0,0];                  % Rates

s_0_2a = [pos_0 ori_0 vel_0 omega_0];

% Force/Moment Calculations
F = ((m*g)/4); %[N]
F_0_2 = F; %[N]
f1 = F; %[N]
f2 = F;  %[N]
f3 = F; %[N]
f4 = F; %[N]
[Lc,Mc,Nc] = controlLMNc(radius,f1,f2,f3,f4,km); % Calculating  control moments

% ODE 45 Function Call
t_span = [0 10]; % Time Span
[t2a,s2a] = ode45(@odefun,t_span,s_0_2a);  

% Plots
plot12(t2a,s2a)
fprintf("Problem 2a: As seen in Figures 1 and 2, the steady hovering flight of our quadrotor is not affected by our aerodynamic forces and moments.\n");
%% Problem (2b)

% Initial Conditions
speed = 5; %[m/s]

Y = -vCoef*speed^2; %[N]
phi = atan(Y/(m*g)); %[rad]

pos_0 = [0,0,-5];                   % Position %[m]
ori_0 = [phi,0,0];                    % Orientation %[rad]
vel_0 = [0,speed*cos(phi),-speed*sin(phi)];% Velocity %[m/s]
omega_0 = [0,0,0];                  % Rates [rad/s]

s_0_2b = [pos_0 ori_0 vel_0 omega_0];

% Force/Moment Calculations
Z_0 = -vCoef * speed * vel_0(3);
F = -m*(g*cos(phi)+Z_0/m);
f1 = F/4; %[N]
f2 = F/4;  %[N]
f3 = F/4; %[N]
f4 = F/4; %[N]
[Lc,Mc,Nc] = controlLMNc(radius,f1,f2,f3,f4,km); % Calculating  control moments

% ODE 45 Function Call
t_span = [0 10]; % Time Span
[t2b,s2b] = ode45(@(t,x) odefun(t,s_0_2b),t_span,s_0_2b);  

% Plots
plot12(t2b,s2b)
fprintf("Problem 2b: The trim thrust values are %d [N]\n",F);
fprintf("Problem 2b: Trim State:");
fprintf(" %f m, %f m, %f m;\n %f rad, %f rad, %f rad;\n %f m/s, %f m/s, %f m/s;\n %f rad/s, %f rad/s, %f rad/s;\n",...
   s2b(1,1), s2b(1,2), s2b(1,3),...
   s2b(1,4), s2b(1,5), s2b(1,6),...
   s2b(1,7), s2b(1,8), s2b(1,9),...
   s2b(1,10), s2b(1,11), s2b(1,12))
%% Problem (2c)

speed = 5; %[m/s]

Y = -vCoef*speed^2;
theta = -atan(Y/(m*g));

pos_0 = [0,0,-5];                   % Position
ori_0 = [0,theta,pi/2];                    % Orientation
vel_0 = [speed*cos(theta),0,speed*sin(theta)];                    % Velocity
omega_0 = [0,0,0];                  % Rates

s_0_2c = [pos_0 ori_0 vel_0 omega_0];

% Force/Moment Calculations
Z_0 = -vCoef * speed * vel_0(3);
F = -m*(g*cos(theta)+Z_0/m);
f1 = F/4; %[N]
f2 = F/4; %[N] 
f3 = F/4; %[N]
f4 = F/4; %[N]
[Lc,Mc,Nc] = controlLMNc(radius,f1,f2,f3,f4,km); % Calculating  control moments

% ODE 45 Function Call
t_span = [0 10]; % Time Span
[t2c,s2c] = ode45(@(t,x) odefun(t,s_0_2c),t_span,s_0_2c);  

% Plots
plot12(t2c,s2c)
fprintf("Problem 2c: The trim thrust values are %d [N]\n",F);
fprintf("Problem 2c: Trim State:");
fprintf(" %f m, %f m, %f m;\n %f rad, %f rad, %f rad;\n %f m/s, %f m/s, %f m/s;\n %f rad/s, %f rad/s, %f rad/s;\n",...
   s2b(1,1), s2b(1,2), s2b(1,3),...
   s2b(1,4), s2b(1,5), s2b(1,6),...
   s2b(1,7), s2b(1,8), s2b(1,9),...
   s2b(1,10), s2b(1,11), s2b(1,12))
%% Problem (3) 
% Load in data file
file = load('RSdata_nocontrol.mat');
Data = file.rt_estim.signals.values;
Time = file.rt_estim.time(:);

% Initial Conditions
pos_0 = [0 0 -5];                  % Position
ori_0 = zeros(1,3);                % Orientation
vel_0 = zeros(1,3);                % Velocity
omega_0 = [0 0 .5];                % Rates
s_0_3 = [pos_0 ori_0 vel_0 omega_0];

% Force/Moment Calculations
Z_0 = -vCoef * speed * vel_0(3);
F =  m*g;
f1 = 0;
f2 = 0; 
f3 = 0;
f4 = 0;
[Lc,Mc,Nc] = controlLMNc(radius,f1,f2,f3,f4,km); % Calculating  control moments

t_span = [0 10];
[t3c,s3c] = ode45(@odefun_level,t_span,s_0_3);  

s3 = [Data(:,1:3),Data(:,4:6),Data(:,7:9),Data(:,10:12)];

% Plots
plot12(Time,s3)
plot12(t3c,s3c)
fprintf("Problem 3: The trim thrust values are %d [N]\n",F);
fprintf("Problem 3: Trim State:");
fprintf(" %f m, %f m, %f m;\n %f rad, %f rad, %f rad;\n %f m/s, %f m/s, %f m/s;\n %f rad/s, %f rad/s, %f rad/s;\n",...
   s3c(1,1), s3c(1,2), s3c(1,3),...
   s3c(1,4), s3c(1,5), s3c(1,6),...
   s3c(1,7), s3c(1,8), s3c(1,9),...
   s3c(1,10), s3c(1,11), s3c(1,12))
fprintf("As can be seen from the physical system in figure 5 compared to the simulated steady hovering flight figure 6,\n")
fprintf("there are some discrepencies between our model and the physical system. From our simulation, we can assume that\n")
fprintf("under ideal simulation conditions, a quad rotor would be stable. Related to problems 1A and 2A, we can see that\n")
fprintf("aerodynamic forces from drag don't have a large effect on stationary craft and it will indeed stay level. However, once moving, our drag forces\n")
fprintf("scale with our velocity to cause more variation in our state variables as shown in figures 3-4. Upon investigating the video, it can be assumed\n")
fprintf("that external forces outside of the scope of our simulation, such as drafts or wind, could be the source of the instability\n")
fprintf("of the quad rotor. But under perfect conditions with no external forces, a quad rotor is stable.\n")


%% Window
fprintf('Figure 1 plots steady level flight at 5m WITHOUT aerodynamics forces.\nFigure 2 plots steady level flight at 5m WITH aerodynamics forces.\n')
fprintf('Figure 3 plots the quadcopter with 5 m/s East velocity with yaw = 0 deg.\nFigure 4 plots the quadcopter with 5 m/s East velocity with yaw = 90 deg.\n')
fprintf('Figure 5 plots the data from the Parrot Mambo drone for stability analysis.\nFigure 6 is the simulation with the initial conditions from the Parrot Mambo Physical data.\n')
%% Functions
function [Lc,Mc,Nc] = controlLMNc(radius,f1,f2,f3,f4,km)
    Lc = (radius/sqrt(2))*(-f1-f2+f3+f4); 
    Mc = (radius/sqrt(2))*(f1-f2-f3+f4); 
    Nc = km*(f1-f2+f3-f4);
end