%TBM Thermal Calculations Code
%Zak Reichenbach
%6/14/2022

%% Linkes
%motor https://www.electricmotorsport.com/me1616-brushless-65hp-liquid-cooled-ipm-motor-24-120v.html

%% House Keeping
clc
clear all
close all

mm2in = 0.0393701; %LINEAR DISTANCE
m2in2 = 1550; %SQUARE DISTANCE
m2in = 39.370;
%% Variables
MotorPower = 5; % [kW]

Qmotor = MotorPower*0.08*1000; %[W] 92% efficiency on motor
Qprop = 0.6*1000;       %for all 6 prop linear actuators based off eff. loss
% Qgear = 0.5*1000;       %safety if gearbox is cooled
QThermal = Qmotor+Qprop ; %+Qgear;
h = 500; % [W/m^2k] https://www.engineersedge.com/heat_transfer/convective_heat_transfer_coefficients__13378.htm 
Tsurf =60+273; %K linear actuator temp 100C for motor
Tinf = 0+273; %estimated water temp at start of loop https://www.thunderstruck-ev.com/me1616.html#:~:text=Coolant%20temp%20leaving%20the%20motor,the%20inlet%20side%20is%20allowable.
%% Calculations

%New goal is to look into convection coeff. with volume flow for the pump
%and the amount of water needed to hold the heat for steady state

%https://www.engineersedge.com/thermodynamics/convective_heat_transfer.htm
D = 0.5; %in
D_m = D/m2in;
r =D/2;

kWat = 0.545; %W/m*k
kCop = 398; %W/m*k

%Qmotor + Qgear + QProp = H*A*(Ts-Tinf);
%This is the area needed for surface area within the reservoir for the
%icebath
A_M = (Qmotor+Qprop)/((Tsurf-Tinf)*(h)); %m^2
A_in = A_M*m2in2;

%Cylinder Surface area 2*PI*r*L 
%r = 1/2D''

Pipelength = A_in/(2*pi*0.25); % Pipe length


%% Dittus-Boelter Equation (Pg. 11 of Convection Heat Transfer Coefficient Estimation by Dr.Harlan H. Bengtson)
%This is to find velocity with a given estimated convection coeff.


%L = V/SA of loop. or diameter of pipe/length of pipe
L = D/Pipelength;
%K = thermal conductivity


% Reynolds Number  = DVp/mu

% MIX SOLUTION https://www.engineeringtoolbox.com/ethylene-glycol-d_146.html

mu = 0.3550*10^-3; %Ns/M^2
rho = 1025; %kg/m^3 ~25% gycol Water mix
Cp = 4.1868*.935*1000; %j/kgk  % ~25% gycol Water mix

%Prandtl Number = Cp*mu/K;

Pr = (Cp) * mu/kWat;

%Nusselt Number = 0.023 * Re^(0.8)*Pr^(0.4)


%hc = (Nu * k)/L Convection Coeff

%LaminarFlow in Circular Tube

Nu = 4.36;


Re = (Nu/(0.023*Pr^0.4))^(1/.8);
%0.6<=PR<=160

Vel = Re*mu/(D_m*rho);

VolFlow = Vel * A_M; % M^3/s

VolFlowPump = VolFlow * 951019; % Gal/hr

%Cp = Specific heat capacity
%K = thermal Conductivity



%Tp =  Temp Plate
%Ta = Temp Air
%Beta =  air thermal expansion coeff
%D = Diameter of Pipe

m = QThermal/(Cp*(10))



