% ASEN 3128 Lab 5 Main
% By Cannon Palmer, Samuel Barrette, Sruthi Bandla, Zak Reichenbach

%clear all;
close all;

%%% aircraft parameter file
recuv_tempest;

tfinal = 100;

TSPAN = [0 tfinal];

wind_inertial = [0;0;0];

%%% set initial conditions
position_inertial0 = [0;0;-600];
euler_angles0 = [0;0;0];
velocity_body0 = [20; 0; 5];
omega_body0 = [0;0;0];

%Initial Conditions for ode45
aircraft_state0 = [position_inertial0; euler_angles0; velocity_body0; omega_body0];

control_input0 = [20*pi/180; 0; 0; .75];

%%% Full sim in ode45
[TOUT,YOUT] = ode45(@(t,y) AircraftEOM(t,y,control_input0,wind_inertial,aircraft_parameters),TSPAN,aircraft_state0,[]);

for i=1:length(TOUT)
    UOUT(i,:) = control_input0';
end

%%% plot results
% PlotAircraftSim(TOUT,YOUT,UOUT,'b')

%% Question 3
v0 = 19;

h0 = 600;

[trim_variables3,fval] = CalculateTrimVariables([v0;h0],aircraft_parameters);

[aircraft_state_trim3,control_surfaces_trim3] = initialTrim([v0;h0],trim_variables3);

[TOUT,YOUT] = ode45(@(t,y) AircraftEOM(t,y,control_surfaces_trim3,wind_inertial,aircraft_parameters),TSPAN,aircraft_state_trim3,[]);

UOUT3 = zeros(length(TOUT),4);

wind = zeros(length(TOUT),3);

for i=1:length(TOUT)
     UOUT3(i,:) = control_surfaces_trim3;
     wind(i,:) = wind_inertial';
end

%PlotAircraftSim(TOUT, YOUT, UOUT3, '--b', wind,0)
 
%% Question 3 Part B
[alpha_trim,de_trim] = CalculateTrimFromStaticStability([v0;h0],aircraft_parameters);

trim_variables3b = [alpha_trim,de_trim,trim_variables3(3)];

[aircraft_state_trim3,control_surfaces_trim3] = initialTrim([v0;h0],trim_variables3b);

[TOUT,YOUT] = ode45(@(t,y) AircraftEOM(t,y,control_surfaces_trim3,wind_inertial,aircraft_parameters),TSPAN,aircraft_state_trim3,[]);

UOUT3 = zeros(length(TOUT),4);

wind = zeros(length(TOUT),3);

for i=1:length(TOUT)
    UOUT3(i,:) = control_surfaces_trim3;
    wind(i,:) = wind_inertial';
end
%PlotAircraftSim(TOUT, YOUT, UOUT3, '-.r', wind,7)

%% Question 4
%Initial Setup
v0 = 21;
h0 = 2000;
tfinal = 150;
TSPAN = [0 tfinal];
wind_inertial = [0;0;0];

%Part a
[trim_variables4,fval] = CalculateTrimVariables([v0;h0],aircraft_parameters);

[aircraft_state_trim4,control_surfaces_trim4] = initialTrim([v0;h0],trim_variables4);

[TOUT,YOUT] = ode45(@(t,y) AircraftEOM(t,y,control_surfaces_trim4,wind_inertial,aircraft_parameters),TSPAN,aircraft_state_trim4,[]);

UOUT3 = zeros(length(TOUT),4);

wind = zeros(length(TOUT),3);

for i=1:length(TOUT)
    UOUT3(i,:) = control_surfaces_trim4;
    wind(i,:) = wind_inertial';
end

PlotAircraftSim(TOUT, YOUT, UOUT3, '--b', wind,14)

%Part b
wind_inertial = [10;10;0];

[TOUT,YOUT] = ode45(@(t,y) AircraftEOM(t,y,control_surfaces_trim4,wind_inertial,aircraft_parameters),TSPAN,aircraft_state_trim4,[]);

UOUT3 = zeros(length(TOUT),4);

wind = zeros(length(TOUT),3);

for i=1:length(TOUT)
    UOUT3(i,:) = control_surfaces_trim4;
    wind(i,:) = wind_inertial';
end

PlotAircraftSim(TOUT, YOUT, UOUT3, '--k', wind,14)

%Part c
wind_inertial = [10;10;0];

Veb = TransformFromInertialToBody(wind_inertial,[0;trim_variables4(1);0]) + TransformFromInertialToBody([v0;0;0],[0;trim_variables4(1);0]);

V = norm(Veb);

alpha0 = trim_variables4(1);

beta0 = asin(Veb(2)/V);

u = V * cos(alpha0) * cos(beta0);

v = V * sin(beta0);

w = V * sin(alpha0) * cos(beta0);

aircraft_state_trim4c = [0,0,-h0,0,alpha0,0,u,v,w,0,0,0]';

[TOUT,YOUT] = ode45(@(t,y) AircraftEOM(t,y,control_surfaces_trim4,wind_inertial,aircraft_parameters),TSPAN,aircraft_state_trim4c,[]);

UOUT3 = zeros(length(TOUT),4);

wind_angles = WindAnglesFromVelocityBody(Veb);

wind = zeros(length(TOUT),3);

for i=1:length(TOUT)
    UOUT3(i,:) = control_surfaces_trim4;
    wind(i,:) = wind_angles;
end

PlotAircraftSim(TOUT, YOUT, UOUT3, '--r', wind,14)

%part d
ca = deg2rad(60); %course angle defined in radians

yaw = atan((10-10*tan(-(pi/2)+ca))/(21*tan(-(pi/2)+ca)));

Veb = TransformFromInertialToBody(wind_inertial,[0;trim_variables4(1);yaw]) + TransformFromInertialToBody([v0;0;0],[0;trim_variables4(1);0]);

V = norm(Veb);

alpha0 = trim_variables4(1);

beta0 = asin(Veb(2)/V);

u = V * cos(alpha0) * cos(beta0);

v = V * sin(beta0);

w = V * sin(alpha0) * cos(beta0);

aircraft_state_trim4d = [0,0,-h0,0,alpha0,yaw,u,v,w,0,0,0]';

[TOUT,YOUT] = ode45(@(t,y) AircraftEOM(t,y,control_surfaces_trim4,wind_inertial,aircraft_parameters),TSPAN,aircraft_state_trim4d,[]);

UOUT3 = zeros(length(TOUT),4);

wind = zeros(length(TOUT),3);

for i=1:length(TOUT)
    UOUT3(i,:) = control_surfaces_trim4;
    wind(i,:) = wind_inertial';
end

PlotAircraftSim(TOUT, YOUT, UOUT3, '--g', wind,14)

%% Plot Details
figure(15)
legend('4.a','4.b','4.c','4.d')

figure(16)
legend('4.a','4.b','4.c','4.d')

figure(17)
legend('4.a','4.b','4.c','4.d')

figure(18)
legend('4.a','4.b','4.c','4.d')

figure(19)
legend('4.a','4.b','4.c','4.d')
grid on

figure(20)
legend('4.a','4.b','4.c','4.d')

figure(21)
legend('4.a','4.b','4.c','4.d')

