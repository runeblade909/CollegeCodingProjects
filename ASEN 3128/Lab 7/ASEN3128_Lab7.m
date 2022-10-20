%% ASEN3128_Lab7.m
clear; close all; clc;

%% Givens
recuv_tempest;
h = 2438.5; %m - Height
V = 22; %m/s - Airspeed
trim_definition = [V, h]; % Trim definition vector

%% Problem 2a
% Find A & B matrices
[trim_variables, fval] = CalculateTrimVariables(trim_definition, aircraft_parameters); % Calculate trim variables
[A_lon, B_lon, A_lat, B_lat] = AircraftLinearModel(trim_definition, trim_variables, aircraft_parameters); % Aircraft linear model

% Root Locus
[numR, denR] = ss2tf(A_lat, B_lat(:,2), [0 0 1 0 0 0], 0);
figure(1);
rlocus(-numR, denR)
hold on
rlocus(numR, denR)
title('Root Locus - Yaw Damper')

%% Problem 2b
% 0 gain Damping Ratio
E_lat = eig(A_lat); %A matrix eigenvalues
xi_lat = sqrt(real(E_lat).^2./(imag(E_lat).^2+real(E_lat).^2)); %Damping ratio
xi_DR_c = xi_lat(4)*1.5; %Increase damping ratio of dutch roll by 50% (target damping ratio)

% From locus plot generated in 2a
k_r = -3.26; %Estimated gain
K = [0 0 k_r 0 0 0]; %Gain vector

%% Problem 2c
% Updated state space matrix parameters
controlSS_DR = A_lat-(B_lat(:,2)*K); %Updated state space matrix
E_DR_control = eig(controlSS_DR); %State space matrix eigenvalues
xi_DR_control = sqrt(real(E_DR_control).^2./(imag(E_DR_control).^2+real(E_DR_control).^2)); %State space matrix damping ratios
omega_n_DR_control = -real(E_DR_control)./xi_DR_control; %Natural frequency

% 0 gain state space matrix parameters
omega_n_lat = -real(E_lat)./xi_lat; %Natural frequency

% Print Results
fprintf('\nDutch Roll Mode:\n                 Uncontrolled: xi = %.4f, omega_n = %.3f.\n                 Controlled:   xi = %.4f, omega_n = %.3f.\n\n', xi_lat(4), omega_n_lat(4), xi_DR_control(4), omega_n_DR_control(4))
fprintf('Roll Mode:\n                 Uncontrolled: omega_n = %.2f.\n                 Controlled:   omega_n = %.2f.\n\n',omega_n_lat(3), omega_n_DR_control(3))
fprintf('Spiral Mode:\n                 Uncontrolled: omega_n = %.2f.\n                 Controlled:   omega_n = %.2f.\n\n',omega_n_lat(6), omega_n_DR_control(6))

%% Problem 2d
% Time
tmax = 10; %s - Maximum simulation time
tspan = [0 tmax]; %s - Time span

% Control surfaces & Initial state
[trim_state, trim_input]= TrimStateAndInput(trim_variables, trim_definition);
x0 = trim_state;
x0(12) = x0(12)+0.015;

% Wind
wind_inertial = [0; 0; 0];

% ODE45
opts = odeset('Events',@phase);
[t, x] = ode45(@(t,x) AircraftEOMControlyaw(t, x, trim_input, wind_inertial, aircraft_parameters, k_r), tspan, x0, opts);

% Plotting
col = 'b';
PlotAircraftSim(t, x, [], wind_inertial, col)

% Uncontrolled ODE45
k_r = 0;
[tu, xu] = ode45(@(tu,xu) AircraftEOMControlyaw(tu, xu, trim_input, wind_inertial, aircraft_parameters, k_r), tspan, x0, opts);

% Plotting
col = 'r';
PlotAircraftSim(tu, xu, [], wind_inertial, col)

%% Phase - ODE options function
function [value, istermainal, direction] = phase(~,x)
%
% Inputs: x = state vector
%         
% Outputs: N/A
%
% Methodology: Determines if the ODE function is terminal. Will stop ODE45
%              if the z value goes above 0.
    
value = x(3);
istermainal = 1;
direction = 1;
end
