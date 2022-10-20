%% ASEN3128_Lab7.m
clear; close all; clc;

%% Givens
recuv_tempest;
h = 2438.5; %m - Height
V = 22; %m/s - Airspeed
trim_definition = [V, h]; % Trim definition vector

%% Problem 4c
% Time
tmax = 100; %s - Maximum simulation time
tspan = [0 tmax]; %s - Time span

[trim_variables, fval] = CalculateTrimVariables(trim_definition, aircraft_parameters); % Calculate trim variables

% Control surfaces & Initial state
[trim_state, trim_input]= TrimStateAndInput(trim_variables, trim_definition);
x0 = trim_state;
x0(12) = x0(12)+0.015;

% Wind
wind_inertial = [0; 0; 0];

%ODE 45 For roll Inner Loop
k_a = -11;
k_p = 1;
opts = odeset('Events',@phase);
[t, x] = ode45(@(tu,xu) AircraftEOMControlRoll(tu, xu, trim_input, wind_inertial, aircraft_parameters, 0, k_a, k_p), tspan, x0, opts);


% % ODE45
% k_r = -70;
% k_a = -100;
% k_p = .1;
% opts = odeset('Events',@phase);
% [t, x] = ode45(@(tu,xu) AircraftEOMControlBoth(tu, xu, trim_input, wind_inertial, aircraft_parameters, 0, k_a,k_p,k_r), tspan, x0, opts);

% Plotting
col = 'b';
PlotAircraftSim(t, x, [], wind_inertial, col)


%ODE 45 For roll Outter Loop
k_a = -11;
k_p = 0;
opts = odeset('Events',@phase);
[tu, xu] = ode45(@(tu,xu) AircraftEOMControlRoll(tu, xu, trim_input, wind_inertial, aircraft_parameters, 0, k_a, k_p), tspan, x0, opts);




% % No yaw ODE45
% k_r = 0;
% k_a = -100;
% k_p = .1;
% [tu, xu] = ode45(@(tu,xu) AircraftEOMControlBoth(tu, xu, trim_input, wind_inertial, aircraft_parameters, 0, k_a,k_p,k_r), tspan, x0, opts);

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
