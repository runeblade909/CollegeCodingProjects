%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% CODE CHALLENGE 7 - Template Script

%

% The purpose of this challenge is to estimate the velocity and kinetic

% energy profile of a falling object. 

%

% To complete the challenge, execute the following steps:

% 1) Set an initial condition velocity

% 2) Set values for constants

% 3) Propagate freefall w/ drag for 20 seconds

% 4) Plot the velocity vs. time

% 5) Calculate the change kinetic energy vs. time

% 6) Plot the change in kinetic energy vs. time

%

% NOTE: DO NOT change any variable names already present in the code.

% 

% Upload your team's script to Gradescope when complete.

% 

% NAME YOUR FILE AS Challenge7_Sec{section number}_Group{group breakout #}.m 

% ***Section numbers are 1 or 2*** 

% EX File Name: Challenge7_Sec1_Group15.m 

%

% STUDENT TEAMMATES

% 1) Austin Coleman

% 2) Andrew Yarbrough

% 3) Austin Sommars

% 4) Zak Reichenbach

% 5) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 

%% Housekeeping

clear variables; close all; clc;

 

%% Set up

m = .3; % [kg]

g = 9.81; % [m/s^2]

rho = 1.225; % [kg/m^3]

Cd = 1.2; % coefficient of drag

A = .0046; % [m^2]

v0 = 0; % [m/s]

 

x0 = 0; % Initial Time - units: [s]

xf = 20; % Final Time - units: [s]

 

 

% Drag = Cd * .5 * (rho * V^2) * A

% F = ma

% F = ?

% ((Cd * .5 * (rho * V^2) * A) - (m*g)) / m = a;

 

 

%% Propagate with ode45

[t,v] = ode45 (@(t,v) g_fun(t,v,m,g,rho,Cd,A), [x0 xf], v0); % Implement ODE45 function.

 

%% Plot Velocity vs. Time

subplot(1,2,1) % Specify location in subplot

plot(t,v);

title('Velocity vs. Time');

xlabel('Time (s)');

ylabel('Velocity (m/s)');

 

 

%% Calculate Kinetic Energy 

KE = (.5 * m) .* (v.^2); % Calculate kinetic energy of plane

 

 

%% Plot Kinetic Energy vs. Time

subplot(1,2,2) % Specify location in subplot

plot(t,KE);

title('Kinetic Energy vs. Time');

xlabel('Time (s)');

ylabel('Kinetic Energy (J)');
