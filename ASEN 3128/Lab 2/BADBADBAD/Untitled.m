%%% 3128 Lab 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Written by: Chesney Boal, Cole MacPherson, Lorenzo Madera
% Group 3: Cole MacPherson, Madisen Frie, Lorenzo Madera 
% ID #s: CB:107519244, CM: 108521329, MF: 107833596, LM: 104733248
% OBJECTIVES:        
% Simulate the EOM for a quadrotor.         
% Investigate the concepts of trim and static stability.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars
close all
clc

%% Problem 1  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create a simulation of a quadrotor: 
% Add control forces and moments from the rotors as discussed in class:

% Initial conditions of quadrotor 
m = 0.068;      % kg
r_dist = 0.060; % m radial dist from CG to prop
k_m = 0.0024;   % N*m/(N) control moment coeff
I_x = 6.8E-5;   % kg*m^2 Body x-axis Moment of Inertia 
I_y = 9.2E-5;   % kg*m^2 Body y-axis Moment of Inertia
I_z = 1.35E-4;  % kg*m^2 Body z-axis Moment of Inertia
nu = 1E-3;      % N/(m/s)^2 Aerodynamic force coefficient 
mu = 2E-6;      % N*m/(rad/s)^2 Aerodynamic moment coefficient
g = 9.81;       % m/s^2 

% set forces in balance 
% THESE ARE THE TRIM THRUSTS
% Initial Control force Z_c0 opposes gravity 
% Z_c0 = mg = 0.6671 = m(w_E_dot - (q*u_E - p*v_E) + g(cos(theta)*cos(phi))
f1 = g/4*m;
f2 = f1;
f3 = f1;
f4 = f1;
constants = [m r_dist k_m I_x I_y I_z nu mu g f1 f2 f3 f4];
disp("Question 1:")
fprintf("The trim thrusts of the rotors for steady hovering flight are: %fN, %fN, %fN, and %fN\n", f1, f2, f3, f4) 
 
% choose a time span 
times = [0.1:100];

% SET INTIAL STATE VECTOR 
% intertial postion
x_data = 1 + zeros(100,1);
y_data = 1 + zeros(100,1);
z_data = -1 + zeros(100,1);
% Euler angles 
phi = zeros(100,1);
theta = zeros(100,1);
psi = zeros(100,1);
% interal vel vectors 
u_E_data = zeros(100,1);
v_E_data = zeros(100,1);
w_E_data = zeros(100,1);
% angular velocity 
p_data = zeros(100,1);
q_data = zeros(100,1);
r_data = zeros(100,1);

% create time span for ODE 
tspan = [times(1,1) times(1,end)]; 

% initial conditions 
x_0 = [x_data(1,1) y_data(1,1) z_data(1,1)...
    phi(1,1) theta(1,1) psi(1,1)...
    u_E_data(1,1) v_E_data(1,1) w_E_data(1,1)...
    p_data(1,1) q_data(1,1) r_data(1,1)];

% run ODE to find trim state 
[t, statevec1] = ode45(@(t,x) objectEOM1(t,x_0,constants), tspan, x_0);
disp("Trim State Vect:")
fprintf("%fkm, %fkm, %fkm;\n %frad, %frad, %frad;\n %fkm/s, %fkm/s, %fkm/s;\n %frad/s, %frad/s, %frad/s;\n",...
   statevec1(1,1), statevec1(1,2), statevec1(1,3),...
   statevec1(1,4), statevec1(1,5), statevec1(1,6),...
   statevec1(1,7), statevec1(1,8), statevec1(1,9),...
   statevec1(1,10), statevec1(1,11), statevec1(1,12))

figure 
sgtitle("Question 1: Including Control F and M, Postion over Time")
subplot(3,1,1)
plot(t,statevec1(:,1))
xlabel("Time [s]")
ylabel("X [m]")
subplot(3,1,2)
plot(t,statevec1(:,2))
xlabel("Time [s]")
ylabel("Y [m]")
subplot(3,1,3)
plot(t,statevec1(:,3))
xlabel("Time [s]")
ylabel("Z [m]")

figure 
sgtitle("Question 1: Including Control F and M, Euler Anlges over Time")
subplot(3,1,1)
plot(t,statevec1(:,4))
xlabel("Time [s]")
ylabel("Roll [rad]")
subplot(3,1,2)
plot(t,statevec1(:,5))
xlabel("Time [s]")
ylabel("Pitch [rad]")
subplot(3,1,3)
plot(t,statevec1(:,6))
xlabel("Time [s]")
ylabel("Yaw [rad]")

figure 
sgtitle("Question 1: Including Control F and M, Velocity over Time in Body Coord")
subplot(3,1,1)
plot(t,statevec1(:,7))
xlabel("Time [s]")
ylabel("uE [m/s]")
subplot(3,1,2)
plot(t,statevec1(:,8))
xlabel("Time [s]")
ylabel("vE [m/s]")
subplot(3,1,3)
plot(t,statevec1(:,9))
xlabel("Time [s]")
ylabel("wE [m/s]")

figure 
sgtitle("Question 1: Including Control F and M, Angular Velocity over Time")
subplot(3,1,1)
plot(t,statevec1(:,10))
xlabel("Time [s]")
ylabel("p [rad/s]")
subplot(3,1,2)
plot(t,statevec1(:,11))
xlabel("Time [s]")
ylabel("q [rad/s]")
subplot(3,1,3)
plot(t,statevec1(:,12))
xlabel("Time [s]")
ylabel("r [rad/s]")

disp("The simulation of Q1 verifies equilibruim motion becuase")
disp("none of the state conditions vary with time.")

%% Problem 2 a 
% Proceed with the simulation 
% Add aerodynamic forces and moments due to drag:

% run the same initial state and constants through the ODE that includes
% aero forces and moments:
[t, statevec2a] = ode45(@(t,x) objectEOM2(t,x_0,constants), tspan, x_0);

disp(" ")
disp("Question 2 part a:")
disp("Trim State Vect:")
fprintf("%fkm, %fkm, %fkm;\n %frad, %frad, %frad;\n %fkm/s, %fkm/s, %fkm/s;\n %frad/s, %frad/s, %frad/s;\n",...
   statevec2a(1,1), statevec2a(1,2), statevec2a(1,3),...
   statevec2a(1,4), statevec2a(1,5), statevec2a(1,6),...
   statevec2a(1,7), statevec2a(1,8), statevec2a(1,9),...
   statevec2a(1,10), statevec2a(1,11), statevec2a(1,12))

figure 
sgtitle("Question 2 a: Including Aero F and M, Postion over Time")
subplot(3,1,1)
plot(t,statevec2a(:,1))
xlabel("Time [s]")
ylabel("X [m]")
subplot(3,1,2)
plot(t,statevec2a(:,2))
xlabel("Time [s]")
ylabel("Y [m]")
subplot(3,1,3)
plot(t,statevec2a(:,3))
xlabel("Time [s]")
ylabel("Z [m]")

figure 
sgtitle("Question 2 a: Including Aero F and M, Euler Anlges over Time")
subplot(3,1,1)
plot(t,statevec2a(:,4))
xlabel("Time [s]")
ylabel("Roll [rad]")
subplot(3,1,2)
plot(t,statevec2a(:,5))
xlabel("Time [s]")
ylabel("Pitch [rad]")
subplot(3,1,3)
plot(t,statevec2a(:,6))
xlabel("Time [s]")
ylabel("Yaw [rad]")

figure 
sgtitle("Question 2 a: Including Aero F and M, Velocity over Time in Body Coord")
subplot(3,1,1)
plot(t,statevec2a(:,7))
xlabel("Time [s]")
ylabel("uE [m/s]")
subplot(3,1,2)
plot(t,statevec2a(:,8))
xlabel("Time [s]")
ylabel("vE [m/s]")
subplot(3,1,3)
plot(t,statevec2a(:,9))
xlabel("Time [s]")
ylabel("wE [m/s]")

figure 
sgtitle("Question 2 a: Including Aero F and M, Angular Velocity over Time")
subplot(3,1,1)
plot(t,statevec2a(:,10))
xlabel("Time [s]")
ylabel("p [rad/s]")
subplot(3,1,2)
plot(t,statevec2a(:,11))
xlabel("Time [s]")
ylabel("q [rad/s]")
subplot(3,1,3)
plot(t,statevec2a(:,12))
xlabel("Time [s]")
ylabel("r [rad/s]")

% statevec 2a [41x12] is all 0 which is the same at statevec1 therefore
% does not alter the trim state 
disp("The simulation of Q2a verifies equilibruim motion becuase")
disp("none of the state conditions vary with time.")


%% Problem 2 b
% Question: Det trim state and rotor thrust trim values
% Knowns:
% constant velocity translation at 5 m/s East
% while maintaining a yaw of 0 deg

% FIRST: use the following equations to solve for phi
X_E = 0;
V_E = 5;
W_E = 0;

% solve for phi_0: trim state roll 
phi_0 = atan((nu * V_E * V_E)/(m*g));
theta_0 = 0;
psi_0 = 0;
u_0 = 0;
v_0 = V_E*cos(phi_0);
w_0 = -V_E*sin(phi_0);
X_0 = 0;
Y_0 = -nu * V_E * v_0;
Z_0 = -nu * V_E * w_0;
L_0_c = 0;
M_0_c = 0;
N_0_c = 0;
Z_0_c = -m*(g*cos(phi_0)+Z_0/m);

x_0_E = 1;
y_0_E = 1;
z_0_E = -1;

x2_data = x_0_E + zeros(100,1); % 0 m
y2_data = y_0_E + zeros(100,1); % 0 m
z2_data = z_0_E + zeros(100,1); % 0 m
% Euler angles 
phi2 = phi_0 + zeros(100,1); % roll = 0
theta2 = theta_0 + zeros(100,1); % pitch = 0
psi2 = psi_0 + zeros(100,1); % yaw = 0
% interal vel vectors 
u_E2_data = u_0 + zeros(100,1); % 0 m/s N
v_E2_data = v_0 + zeros(100,1); % 5 m/s E
w_E2_data = w_0 + zeros(100,1); % 0 m/s down
% angular velocity 
p_2data = p_data; % 0 rad/s 
q_2data = q_data; % 0 rad/s 
r_2data = r_data; % 0 rad/s 

% change forces to balance 
f1 = -Z_0_c/4;
f2 = -Z_0_c/4;
f3 = -Z_0_c/4;
f4 = -Z_0_c/4;
constants2 = constants;
constants2(1,10:13) = [f1 f2 f3 f4];

disp(" ")
disp("Question 2 part b:")
fprintf("The trim thrusts of the rotors for 5 m/s E & yaw = 0 deg: %fN, %fN, %fN, and %fN\n", f1, f2, f3, f4) 
 
% create time span for ODE 
tspan = [times(1,1) times(1,end)]; 

x_02 = [x2_data(1,1) y2_data(1,1) z2_data(1,1)...
    phi2(1,1) theta2(1,1) psi2(1,1)...
    u_E2_data(1,1) v_E2_data(1,1) w_E2_data(1,1)...
    p_2data(1,1) q_2data(1,1) r_2data(1,1)];

[t, statevec2b] = ode45(@(t,x) objectEOM2(t,x_02,constants2), tspan, x_02);

disp("Trim State Vect:")
fprintf("%fkm, %fkm, %fkm;\n %frad, %frad, %frad;\n %fkm/s, %fkm/s, %fkm/s;\n %frad/s, %frad/s, %frad/s;\n",...
   statevec2b(1,1), statevec2b(1,2), statevec2b(1,3),...
   statevec2b(1,4), statevec2b(1,5), statevec2b(1,6),...
   statevec2b(1,7), statevec2b(1,8), statevec2b(1,9),...
   statevec2b(1,10), statevec2b(1,11), statevec2b(1,12))

figure 
sgtitle("Question 2 b: 5 m/s East, yaw of 0 deg, Postion over Time")
subplot(3,1,1)
plot(t,statevec2b(:,1))
xlabel("Time [s]")
ylabel("X [m]")
subplot(3,1,2)
plot(t,statevec2b(:,2))
xlabel("Time [s]")
ylabel("Y [m]")
subplot(3,1,3)
plot(t,statevec2b(:,3))
xlabel("Time [s]")
ylabel("Z [m]")
ylim([-1.5 -0.5])

figure 
sgtitle("Question 2 b: 5 m/s East, yaw of 0 deg, Euler Anlges over Time")
subplot(3,1,1)
plot(t,statevec2b(:,4))
xlabel("Time [s]")
ylabel("Roll [rad]")
ylim([0.00 0.05])
subplot(3,1,2)
plot(t,statevec2b(:,5))
xlabel("Time [s]")
ylabel("Pitch [rad]")
subplot(3,1,3)
plot(t,statevec2b(:,6))
xlabel("Time [s]")
ylabel("Yaw [rad]")

figure 
sgtitle("Question 2 b: 5 m/s East, yaw of 0 deg, Velocity over Time in Body Coord")
subplot(3,1,1)
plot(t,statevec2b(:,7))
xlabel("Time [s]")
ylabel("uE [m/s]")
subplot(3,1,2)
plot(t,statevec2b(:,8))
xlabel("Time [s]")
ylabel("vE [m/s]")
ylim([4.99 5.00])
subplot(3,1,3)
plot(t,statevec2b(:,9))
xlabel("Time [s]")
ylabel("wE [m/s]")
ylim([-0.5 0])

figure 
sgtitle("Question 2 b: 5 m/s East, yaw of 0 deg, Angular Velocity over Time")
subplot(3,1,1)
plot(t,statevec2b(:,10))
xlabel("Time [s]")
ylabel("p [rad/s]")
subplot(3,1,2)
plot(t,statevec2b(:,11))
xlabel("Time [s]")
ylabel("q [rad/s]")
subplot(3,1,3)
plot(t,statevec2b(:,12))
xlabel("Time [s]")
ylabel("r [rad/s]")

disp("The simulation of Q2b verifies equilibruim motion becuase")
disp("none of the state conditions vary with time.")

%% Problem 2 c
% Question: What changes in the trim state if:
% Known:
% yaw of 90 deg is to be maintained 
% translating 5 m/s East

X_Ec = 0;
V_Ec = 5;
W_Ec = 0;
phi_0c = 0;
theta_0c = -atan((nu * V_Ec * V_Ec)/(m*g));
psi_0c = 90*pi/180;
u_0c = V_Ec*cos(theta_0c);
v_0c = 0;
w_0c = V_Ec*sin(theta_0c);
X_0c = -nu * V_Ec * u_0c;
Y_0c = 0;
Z_0c = -nu * V_Ec * w_0c;
L_0_c_c = 0;
M_0_c_c = 0;
N_0_c_c = 0;
Z_0_c_c = -m*(g*cos(theta_0c) + Z_0c/m);

x_0_Ec = 1;
y_0_Ec = 1;
z_0_Ec = -1;

x3_data = x_0_Ec + zeros(100,1); % 0 m
y3_data = y_0_Ec + zeros(100,1); % 0 m
z3_data = z_0_Ec + zeros(100,1); % 0 m
% Euler angles 
phi3 = phi_0c + zeros(100,1); % roll = 0
theta3 = theta_0c + zeros(100,1); % pitch = 0
psi3 = psi_0c + zeros(100,1); % yaw = 0
% interal vel vectors 
u_E3_data = u_0c + zeros(100,1); % 0 m/s N
v_E3_data = v_0c + zeros(100,1); % 5 m/s E
w_E3_data = w_0c + zeros(100,1); % 0 m/s down
% angular velocity 
p_3data = p_data; % 0 rad/s 
q_3data = q_data; % 0 rad/s 
r_3data = r_data; % 0 rad/s 

% change forces to balance 
f1 = -Z_0_c_c/4;
f2 = -Z_0_c_c/4;
f3 = -Z_0_c_c/4;
f4 = -Z_0_c_c/4;
constants3 = constants;
constants3(1,10:13) = [f1 f2 f3 f4];

disp(" ")
disp("Question 2 part c:")
fprintf("The trim thrusts of the rotors for 5 m/s E & yaw = 90 deg: %fN, %fN, %fN, and %fN\n", f1, f2, f3, f4) 

% create time span for ODE 
tspan = [times(1,1) times(1,end)]; 

x_03 = [x3_data(1,1) y3_data(1,1) z3_data(1,1)...
    phi3(1,1) theta3(1,1) psi3(1,1)...
    u_E3_data(1,1) v_E3_data(1,1) w_E3_data(1,1)...
    p_3data(1,1) q_3data(1,1) r_3data(1,1)];

[t, statevec2c] = ode45(@(t,x) objectEOM2(t,x_03,constants3), tspan, x_03);

disp("Trim State Vect:")
fprintf("%fkm, %fkm, %fkm;\n %frad, %frad, %frad;\n %fkm/s, %fkm/s, %fkm/s;\n %frad/s, %frad/s, %frad/s;\n",...
   statevec2c(1,1), statevec2c(1,2), statevec2c(1,3),...
   statevec2c(1,4), statevec2c(1,5), statevec2c(1,6),...
   statevec2c(1,7), statevec2c(1,8), statevec2c(1,9),...
   statevec2c(1,10), statevec2c(1,11), statevec2c(1,12))

disp(" ")
disp("The trim state Euler angles and velocity change.")
disp("The trim state roll of 2b has now become the negative pitch angle in 2c.")
disp("The trim state also now requires the 90 degree yaw in 2c.")
disp("The Y component of the trim state velocity of 2b is now the X component in 2c.")

figure 
sgtitle("Question 2 c: 5 m/s East, yaw of 90 deg, Postion over Time")
subplot(3,1,1)
plot(t,statevec2c(:,1))
xlabel("Time [s]")
ylabel("X [m]")
ylim([0.5 1.5])
subplot(3,1,2)
plot(t,statevec2c(:,2))
xlabel("Time [s]")
ylabel("Y [m]")
ylim([0 500])
subplot(3,1,3)
plot(t,statevec2c(:,3))
xlabel("Time [s]")
ylabel("Z [m]")
ylim([-1.5 -0.5])

figure 
sgtitle("Question 2 c: 5 m/s East, yaw of 90 deg, Euler Anlges over Time")
subplot(3,1,1)
plot(t,statevec2c(:,4))
xlabel("Time [s]")
ylabel("Roll [rad]")
subplot(3,1,2)
plot(t,statevec2c(:,5))
xlabel("Time [s]")
ylabel("Pitch [rad]")
ylim([-0.05 0.00])
subplot(3,1,3)
plot(t,statevec2c(:,6))
xlabel("Time [s]")
ylabel("Yaw [rad]")
ylim([1 2])

figure 
sgtitle("Question 2 c: 5 m/s East, yaw of 90 deg, Velocity over Time in Body Coord")
subplot(3,1,1)
plot(t,statevec2c(:,7))
xlabel("Time [s]")
ylabel("uE [m/s]")
ylim([4.99 5.00])
subplot(3,1,2)
plot(t,statevec2c(:,8))
xlabel("Time [s]")
ylabel("vE [m/s]")
subplot(3,1,3)
plot(t,statevec2c(:,9))
xlabel("Time [s]")
ylabel("wE [m/s]")
ylim([-0.5 0])

figure 
sgtitle("Question 2 c: 5 m/s East, yaw of 90 deg, Angular Velocity over Time")
subplot(3,1,1)
plot(t,statevec2c(:,10))
xlabel("Time [s]")
ylabel("p [rad/s]")
subplot(3,1,2)
plot(t,statevec2c(:,11))
xlabel("Time [s]")
ylabel("q [rad/s]")
subplot(3,1,3)
plot(t,statevec2c(:,12))
xlabel("Time [s]")
ylabel("r [rad/s]")

disp("The simulation of Q2c verifies equilibruim motion becuase")
disp("none of the state conditions vary with time.")

%% Problem 3 

% SET INTIAL STATE VECTOR 
% this intial state vector is equal to the same state as Q1 or 2a 
% intertial postion
x_data3a = 1 + zeros(100,1);
y_data3a = 1 + zeros(100,1);
z_data3a = -1 + zeros(100,1);
% Euler angles 
phi3a = zeros(100,1);
theta3a = zeros(100,1);
psi3a = zeros(100,1);
% interal vel vectors 
u_E_data3a = zeros(100,1);
v_E_data3a = zeros(100,1);
w_E_data3a = zeros(100,1);
% angular velocity 
p_data3a = zeros(100,1);
q_data3a = zeros(100,1);
r_data3a = zeros(100,1);

[t, statevec3] = ode45(@(t,x) objectEOM3(t,x_0,constants), tspan, x_0);

disp(" ")
disp("Question 3:")
disp("Trim State Vect without Control Moments or Forces:")
fprintf("%fkm, %fkm, %fkm;\n %frad, %frad, %frad;\n %fkm/s, %fkm/s, %fkm/s;\n %frad/s, %frad/s, %frad/s;\n",...
   statevec3(1,1), statevec3(1,2), statevec3(1,3),...
   statevec3(1,4), statevec3(1,5), statevec3(1,6),...
   statevec3(1,7), statevec3(1,8), statevec3(1,9),...
   statevec3(1,10), statevec3(1,11), statevec3(1,12))

figure 
sgtitle("Question 3: No Control F and M, Postion over Time")
subplot(3,1,1)
plot(t,statevec3(:,1))
xlabel("Time [s]")
ylabel("X [m]")
subplot(3,1,2)
plot(t,statevec3(:,2))
xlabel("Time [s]")
ylabel("Y [m]")
subplot(3,1,3)
plot(t,statevec3(:,3))
xlabel("Time [s]")
ylabel("Z [m]")

figure 
sgtitle("Question 3: No Control F and M, Euler Anlges over Time")
subplot(3,1,1)
plot(t,statevec3(:,4))
xlabel("Time [s]")
ylabel("Roll [rad]")
subplot(3,1,2)
plot(t,statevec3(:,5))
xlabel("Time [s]")
ylabel("Pitch [rad]")
subplot(3,1,3)
plot(t,statevec3(:,6))
xlabel("Time [s]")
ylabel("Yaw [rad]")

figure 
sgtitle("Question 3: No Control F and M, Velocity over Time in Body Coord")
subplot(3,1,1)
plot(t,statevec3(:,7))
xlabel("Time [s]")
ylabel("uE [m/s]")
subplot(3,1,2)
plot(t,statevec3(:,8))
xlabel("Time [s]")
ylabel("vE [m/s]")
subplot(3,1,3)
plot(t,statevec3(:,9))
xlabel("Time [s]")
ylabel("wE [m/s]")

figure 
sgtitle("Question 3: No Control F and M, Angular Velocity over Time")
subplot(3,1,1)
plot(t,statevec3(:,10))
xlabel("Time [s]")
ylabel("p [rad/s]")
subplot(3,1,2)
plot(t,statevec3(:,11))
xlabel("Time [s]")
ylabel("q [rad/s]")
subplot(3,1,3)
plot(t,statevec3(:,12))
xlabel("Time [s]")
ylabel("r [rad/s]")

disp("  ")
disp("Through the behavior of the hardware demonstration video,")
disp("the quadrotor may seem to be in equilibrum, but it is not stable.")
disp("The quadrotor takes an immediate dive near the end of flight.")
disp("When hovering, though, the quadrotor is in equilibrum.")
disp("Through the simulation of Question 3 plots, the system is in equilibrium.")
disp("This simulation does not include disturbances to determine that the flight is stable.")

%% TABLE OF PARTICIPATION!!! INSERT HERE

%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%


%% ODE call prob 1: Includes Control Forces and Moments, but not Aerodynamic F and Ms
function [state_dot] = objectEOM1(t,x,constants)
% Function for ode45 that models the motion of the quadcopter given the
% intial state vector of teh quadcopter and a vector of constants that
% defines the characteristics of the quadcopters. Does not include
% aerodynamic forces and moments
%
% Inputs:  t         - time
%          x         - state vector
%          constants - vector of constants that define the characteristics
%                       of the quadcopter
%
% Outputs: state_dot - the derivative of the state vector

% take state variables out
x_pos = x(:,1);
y = x(:,2);
z = x(:,3);
phi = x(:,4);
theta = x(:,5);
psi = x(:,6);
u = x(:,7);
v = x(:,8);
w = x(:,9);
p = x(:,10);
q = x(:,11);
r = x(:,12);

% take constants out 
m = constants(1,1);
r_dist = constants(1,2);
k_m = constants(1,3);
Ix = constants(1,4);
Iy = constants(1,5);
Iz = constants(1,6);
nu = constants(1,7);
mu = constants(1,8);
g = constants(1,9);
f1 = constants(1,10);
f2 = constants(1,11);
f3 = constants(1,12);
f4 = constants(1,13);

% 1-3
dxdt = [cos(theta)*cos(psi) sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);...
    cos(theta)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);...
    -sin(theta) sin(phi)*cos(theta) cos(phi)*cos(theta);];

x_E_dot = dxdt*x(:,7:9)';

% 4-6
dEulerdt = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);...
        0 cos(phi) -sin(phi);...
        0 sin(phi)*sec(theta) cos(phi)*sec(theta)];

euler_dot = dEulerdt*x(:,10:12)';

% 7-9
inert_vel1 = [r*v - q*w; p*w - r*u; q*u - p*v];
inert_vel2 = [-sin(theta); cos(theta)*sin(phi); cos(theta)*cos(phi)];
f_ctrl = [0; 0; -f1-f2-f3-f4];

inert_vel_dot = inert_vel1 + g*inert_vel2 + 1/m*f_ctrl; 
m_ctrl = [r_dist/sqrt(2)*(-f1-f2+f3+f4); r_dist/sqrt(2)*(f1-f2-f3+f4); k_m*(f1-f2+f3-f4)];
ang_vel1 = [(Iy-Iz)/Ix*q*r; (Iz-Ix)/Iy*p*r; (Ix-Iy)/Iz*p*q];
ang_vel3 = [m_ctrl(1,1)/Ix; m_ctrl(2,1)/Iy; m_ctrl(3,1)/Iz];

ang_vel_dot = ang_vel1 + ang_vel3; %+ ang_vel2

state_dot = [x_E_dot' euler_dot' inert_vel_dot' ang_vel_dot']';

end 

%% ODE call prob 2: Includes Control Forces and Moments AND Aerodynamic F and Ms
function [state_dot] = objectEOM2(t,x,constants)
% Function for ode45 that models the motion of the quadcopter given the
% intial state vector of teh quadcopter and a vector of constants that
% defines the characteristics of the quadcopters.
%
% Inputs:  t         - time
%          x         - state vector
%          constants - vector of constants that define the characteristics
%                       of the quadcopter
%
% Outputs: state_dot - the derivative of the state vector

x_pos = x(:,1);
y = x(:,2);
z = x(:,3);
phi = x(:,4);
theta = x(:,5);
psi = x(:,6);
u = x(:,7);
v = x(:,8);
w = x(:,9);
p = x(:,10);
q = x(:,11);
r = x(:,12);

m = constants(1,1);
r_dist = constants(1,2);
k_m = constants(1,3);
Ix = constants(1,4);
Iy = constants(1,5);
Iz = constants(1,6);
nu = constants(1,7);
mu = constants(1,8);
g = constants(1,9);
f1 = constants(1,10);
f2 = constants(1,11);
f3 = constants(1,12);
f4 = constants(1,13);

% 1-3
dxdt = [cos(theta)*cos(psi) sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);...
    cos(theta)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);...
    -sin(theta) sin(phi)*cos(theta) cos(phi)*cos(theta);];

x_E_dot = dxdt*x(:,7:9)';

% 4-6
dEulerdt = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);...
        0 cos(phi) -sin(phi);...
        0 sin(phi)*sec(theta) cos(phi)*sec(theta)];

euler_dot = dEulerdt*x(:,10:12)';

% 7-9
inert_vel1 = [r*v - q*w; p*w - r*u; q*u - p*v];
inert_vel2 = [-sin(theta); cos(theta)*sin(phi); cos(theta)*cos(phi)];
% f_ctrl [0 0 Zc]]'
f_ctrl = [0; 0; -f1-f2-f3-f4];
f_aero = -nu*sqrt(u^2+v^2+w^2)*x(:,7:9)';

inert_vel_dot = inert_vel1 + g*inert_vel2 + 1/m*f_aero + 1/m*f_ctrl;

% 10-12
m_aero = -mu*sqrt(p^2+q^2+r^2)*x(:,10:12)'; % [L M N]'
% m_ctrl = [Lc Mc Nc]'
m_ctrl = [r_dist/sqrt(2)*(-f1-f2+f3+f4); r_dist/sqrt(2)*(f1-f2-f3+f4); k_m*(f1-f2+f3-f4)];
ang_vel1 = [(Iy-Iz)/Ix*q*r; (Iz-Ix)/Iy*p*r; (Ix-Iy)/Iz*p*q];
ang_vel2 = [m_aero(1,1)/Ix; m_aero(2,1)/Iy; m_aero(3,1)/Iz];
ang_vel3 = [m_ctrl(1,1)/Ix; m_ctrl(2,1)/Iy; m_ctrl(3,1)/Iz];

ang_vel_dot = ang_vel1 + ang_vel2 + ang_vel3;

state_dot = [x_E_dot' euler_dot' inert_vel_dot' ang_vel_dot']'

end 

%% ODE call prob 3: Includes Aerodynamic Forces and Moments, but not Control F and Ms
function [state_dot] = objectEOM3(t,x,constants)
% Function for ode45 that models the motion of the quadcopter given the
% intial state vector of teh quadcopter and a vector of constants that
% defines the characteristics of the quadcopters. Does not include control
% moments and forces
%
% Inputs:  t         - time
%          x         - state vector
%          constants - vector of constants that define the characteristics
%                       of the quadcopter
%
% Outputs: state_dot - the derivative of the state vector

x_pos = x(:,1);
y = x(:,2);
z = x(:,3);
phi = x(:,4);
theta = x(:,5);
psi = x(:,6);
u = x(:,7);
v = x(:,8);
w = x(:,9);
p = x(:,10);
q = x(:,11);
r = x(:,12);

m = constants(1,1);
r_dist = constants(1,2);
k_m = constants(1,3);
Ix = constants(1,4);
Iy = constants(1,5);
Iz = constants(1,6);
nu = constants(1,7);
mu = constants(1,8);
g = constants(1,9);
f1 = constants(1,10);
f2 = constants(1,11);
f3 = constants(1,12);
f4 = constants(1,13);

% 1-3
dxdt = [cos(theta)*cos(psi) sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);...
    cos(theta)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);...
    -sin(theta) sin(phi)*cos(theta) cos(phi)*cos(theta);];

x_E_dot = dxdt*x(:,7:9)';

% 4-6
dEulerdt = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);...
        0 cos(phi) -sin(phi);...
        0 sin(phi)*sec(theta) cos(phi)*sec(theta)];

euler_dot = dEulerdt*x(:,10:12)';

% 7-9
inert_vel1 = [r*v - q*w; p*w - r*u; q*u - p*v];
inert_vel2 = [-sin(theta); cos(theta)*sin(phi); cos(theta)*cos(phi)];
% f_ctrl [0 0 Zc]]'
%f_ctrl = [0; 0; -f1-f2-f3-f4];
f_aero = -nu*sqrt(u^2+v^2+w^2)*x(:,7:9)';

inert_vel_dot = inert_vel1 + g*inert_vel2 + 1/m*f_aero; % + 1/m*f_ctrl

% 10-12
m_aero = -mu*sqrt(p^2+q^2+r^2)*x(:,10:12)'; % [L M N]'
% m_ctrl = [Lc Mc Nc]'
% m_ctrl = [r_dist/sqrt(2)*(-f1-f2+f3+f4); r_dist/sqrt(2)*(f1-f2-f3+f4); k_m*(f1-f2+f3-f4)];
ang_vel1 = [(Iy-Iz)/Ix*q*r; (Iz-Ix)/Iy*p*r; (Ix-Iy)/Iz*p*q];
ang_vel2 = [m_aero(1,1)/Ix; m_aero(2,1)/Iy; m_aero(3,1)/Iz];
% ang_vel3 = [m_ctrl(1,1)/Ix; m_ctrl(2,1)/Iy; m_ctrl(3,1)/Iz];

ang_vel_dot = ang_vel1 + ang_vel2; % + ang_vel3;

state_dot = [x_E_dot' euler_dot' inert_vel_dot' ang_vel_dot']';

end 