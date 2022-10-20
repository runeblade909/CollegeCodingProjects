%% ASEN 3128 - LAb 3 - Main
% Script to compare the linearized and non-linearized models of a quadrotor
% in steady equilibrium flight. Additionally, control moments and forces
% are to be added to model how a quadrotor remains in steady hover flight
%
% Author: Cole MacPherson
% Collaborators: S. Packard, D. Wolfe, D. Zhao
% Date: 5th Mar 2021

%% Housekeeping
clc;
clear;
close all;
tic

%% declare constants
m = 0.068; % mass of the quadrotor [kg]
R = 0.06; % radial distance from CG to propeller [m]
k_m = 0.0024; % control moment coefficient [N*m/N]
I_x = 6.8e-5; % x-axis moment of inertia [kg*m^2]
I_y = 9.2e-5; % y-axis moment of inertia [kg*m^2]
I_z = 1.35e-4; % z-axis moment of inertia [kg*m^2]
nu = 1e-3; % aerodynamic force coefficient [N/(m/s)^2]
mu = 2e-6; % aerodynamic moment coefficient [N*m/(rad/s)^2]
g = 9.81; % graviational constant [m/s^2]

Z_c = -m*g; % control force
L_c = 0; % x control moment
M_c = 0; % y contorl moment
N_c = 0; % z control moment

tspan = [0 5]; % time to integrate over

%% initialize figure information
col = 'b';
fig_a = [1,2,3,4,5,6]; % figure a number vector
fig_b = [1,2,3,4,5,6] + 6; % figure b number vector
fig_c = [1,2,3,4,5,6] + 2*6; % figure c number vector
fig_d = [1,2,3,4,5,6] + 3*6; % figure d number vector
fig_e = [1,2,3,4,5,6] + 4*6; % figure e number vector
fig_f = [1,2,3,4,5,6] + 5*6; % figure f number vector

%% 3a
state_vec_0 = [0 0 0 deg2rad(5) 0 0 0 0 0 0 0 0]';
[t_a,state_vec_a] = ode45(@(t_a,state_vec_a) quadrotorODE(t_a,state_vec_a,m,I_x,I_y,I_z,nu,mu,Z_c,L_c,M_c,N_c),tspan,state_vec_0);

PlotAircraftSim(t_a,state_vec_a,ones(length(t_a),4)*Z_c/4,fig_a,col,'3.a')

%% 3b
state_vec_0 = [0 0 0 0 deg2rad(5) 0 0 0 0 0 0 0]';
[t_b,state_vec_b] = ode45(@(t_b,state_vec_b) quadrotorODE(t_b,state_vec_b,m,I_x,I_y,I_z,nu,mu,Z_c,L_c,M_c,N_c),tspan,state_vec_0);

PlotAircraftSim(t_b,state_vec_b,ones(length(t_b),4)*Z_c/4,fig_b,col,'3.b')

%% 3c
state_vec_0 = [0 0 0 0 0 deg2rad(5) 0 0 0 0 0 0]';
[t_c,state_vec_c] = ode45(@(t_c,state_vec_c) quadrotorODE(t_c,state_vec_c,m,I_x,I_y,I_z,nu,mu,Z_c,L_c,M_c,N_c),tspan,state_vec_0);

PlotAircraftSim(t_c,state_vec_c,ones(length(t_c),4)*Z_c/4,fig_c,col,'3.c')

%% 3d
[p,q,r] = rollrate2pqr(0,0,0,0.1,0,0);
state_vec_0 = [0 0 0 0 0 0 0 0 0 p q r]';
[t_d,state_vec_d] = ode45(@(t_d,state_vec_d) quadrotorODE(t_d,state_vec_d,m,I_x,I_y,I_z,nu,mu,Z_c,L_c,M_c,N_c),tspan,state_vec_0);

PlotAircraftSim(t_d,state_vec_d,ones(length(t_d),4)*Z_c/4,fig_d,col,'3.d')

%% 3e
[p,q,r] = rollrate2pqr(0,0,0,0,0.1,0);
state_vec_0 = [0 0 0 0 0 0 0 0 0 p q r]';
[t_e,state_vec_e] = ode45(@(t_e,state_vec_e) quadrotorODE(t_e,state_vec_e,m,I_x,I_y,I_z,nu,mu,Z_c,L_c,M_c,N_c),tspan,state_vec_0);

PlotAircraftSim(t_e,state_vec_e,ones(length(t_e),4)*Z_c/4,fig_e,col,'3.e')

%% 3f
[p,q,r] = rollrate2pqr(0,0,0,0,0,0.1);
state_vec_0 = [0 0 0 0 0 0 0 0 0 p q r]';
[t_f,state_vec_f] = ode45(@(t_f,state_vec_f) quadrotorODE(t_f,state_vec_f,m,I_x,I_y,I_z,nu,mu,Z_c,L_c,M_c,N_c),tspan,state_vec_0);

PlotAircraftSim(t_f,state_vec_f,ones(length(t_f),4)*Z_c/4,fig_f,col,'3.f')

%% Redefine figure line color for problem 4
col = 'r';

%% 4a
t_final = 5;
state_vec_0 = [0 0 0 deg2rad(5) 0 0 0 0 0 0 0 0]';
[state_vec_a,t_a] = linearizedEOM(state_vec_0,t_final);

PlotAircraftSim(t_a,state_vec_a,ones(length(t_a),4)*Z_c/4,fig_a,col,'3.a & 4.a')

%% 4b
state_vec_0 = [0 0 0 0 deg2rad(5) 0 0 0 0 0 0 0]';
[state_vec_b,t_b] = linearizedEOM(state_vec_0,t_final);

PlotAircraftSim(t_b,state_vec_b,ones(length(t_b),4)*Z_c/4,fig_b,col,'3.b & 4.b')

%% 4c
state_vec_0 = [0 0 0 0 0 deg2rad(5) 0 0 0 0 0 0]';
[state_vec_c,t_c] = linearizedEOM(state_vec_0,t_final);

PlotAircraftSim(t_c,state_vec_c,ones(length(t_c),4)*Z_c/4,fig_c,col,'3.c & 4.c')

%% 4d
[p,q,r] = rollrate2pqr(0,0,0,0.1,0,0);
state_vec_0 = [0 0 0 0 0 0 0 0 0 p q r]';
[state_vec_d,t_d] = linearizedEOM(state_vec_0,t_final);

PlotAircraftSim(t_d,state_vec_d,ones(length(t_d),4)*Z_c/4,fig_d,col,'3.d & 4.d')

%% 4e
[p,q,r] = rollrate2pqr(0,0,0,0,0.1,0);
state_vec_0 = [0 0 0 0 0 0 0 0 0 p q r]';
[state_vec_e,t_e] = linearizedEOM(state_vec_0,t_final);

PlotAircraftSim(t_e,state_vec_e,ones(length(t_e),4)*Z_c/4,fig_e,col,'3.e & 4.e')

%% 4f
[p,q,r] = rollrate2pqr(0,0,0,0,0,0.1);
state_vec_0 = [0 0 0 0 0 0 0 0 0 p q r]';
[state_vec_f,t_f] = linearizedEOM(state_vec_0,t_final);

PlotAircraftSim(t_f,state_vec_f,ones(length(t_f),4)*Z_c/4,fig_f,col,'3.f & 4.f')

%% 5d
[p,q,r] = rollrate2pqr(0,0,0,0.1,0,0);
state_vec_0 = [0 0 0 0 0 0 0 0 0 p q r Z_c -0.004*p -0.004*q -0.004*r]';
[t_d,state_vec_d] = ode45(@(t_d,state_vec_d) quadrotorODE_controlled(t_d,state_vec_d,m,I_x,I_y,I_z,nu,mu,Z_c,L_c,M_c,N_c),tspan,state_vec_0);

% compute forces
F_5d = ComputeMotorForces(state_vec_d(:,13),state_vec_d(:,14),state_vec_d(:,15),state_vec_d(:,16),R,k_m);

% plot force
figure
subplot(3,1,1)
plot(t_d,F_5d(1,:),'linewidth',2); hold on;
plot(t_d,F_5d(2,:),'linewidth',2); hold on;
plot(t_d,F_5d(3,:),'linewidth',2); hold on;
plot(t_d,F_5d(4,:),'linewidth',2); hold on;
title(['p = ' num2str(p)]);
legend('motor 1','motor 2','motor 3','motor 4');
hold on;

%% 5e
[p,q,r] = rollrate2pqr(0,0,0,0,0.1,0);
state_vec_0 = [0 0 0 0 0 0 0 0 0 p q r Z_c -0.004*p -0.004*q -0.004*r]';
[t_e,state_vec_e] = ode45(@(t_e,state_vec_e) quadrotorODE_controlled(t_e,state_vec_e,m,I_x,I_y,I_z,nu,mu,Z_c,L_c,M_c,N_c),tspan,state_vec_0);

% compute forces
F_5e = ComputeMotorForces(state_vec_e(:,13),state_vec_e(:,14),state_vec_e(:,15),state_vec_e(:,16),R,k_m);

% plot force
subplot(3,1,2)
plot(t_e,F_5e(1,:),'linewidth',2); hold on;
plot(t_e,F_5e(2,:),'linewidth',2); hold on;
plot(t_e,F_5e(3,:),'linewidth',2); hold on;
plot(t_e,F_5e(4,:),'linewidth',2); hold on;
title(['q = ' num2str(q)]);
legend('motor 1','motor 2','motor 3','motor 4');
hold on;

%% 5f
[p,q,r] = rollrate2pqr(0,0,0,0,0,0.1);
state_vec_0 = [0 0 0 0 0 0 0 0 0 p q r Z_c -0.004*p -0.004*q -0.004*r]';
[t_f,state_vec_f] = ode45(@(t_f,state_vec_f) quadrotorODE_controlled(t_f,state_vec_f,m,I_x,I_y,I_z,nu,mu,Z_c,L_c,M_c,N_c),tspan,state_vec_0);

% compute forces
F_5f = ComputeMotorForces(state_vec_f(:,13),state_vec_f(:,14),state_vec_f(:,15),state_vec_f(:,16),R,k_m);

% plot force
subplot(3,1,3)
plot(t_f,F_5f(1,:),'linewidth',2); hold on;
plot(t_f,F_5f(2,:),'linewidth',2); hold on;
plot(t_f,F_5f(3,:),'linewidth',2); hold on;
plot(t_f,F_5f(4,:),'linewidth',2); hold on;
title(['r = ' num2str(r)]);
legend('motor 1','motor 2','motor 3','motor 4');
hold on;

%% Plot Housekeeping
% make legends for each plot
for N = 1:6
    figure(6*N-5)
    legend('Non-Linear Path','Linearized Path','location','best');
    legend('boxoff');
    
    figure(6*N-4)
    legend('Non-Linear Path','Linearized Path','location','best');
    legend('boxoff');
    
    figure(6*N-3)
    legend('Non-Linear Path','Linearized Path','location','best');
    legend('boxoff');
    
    figure(6*N-2)
    legend('Non-Linear Path','Linearized Path','location','best');
    legend('boxoff');
    
    figure(6*N-1)
    legend('Non-Linear Path','Linearized Path','location','best');
    legend('boxoff');
    
    figure(6*N)
    legend('Non-Linear Path','Start','End','Linearized Path');
end

%% End Housekeeping