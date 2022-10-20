%Zak Reichenbach
%2003 Lab 2 Individual Part

clc
clear
close all

%Intake appendix 2 data (Column 1 AoA, Column 2 Cl)
load('app2data.mat')
Cl = app2data(:,2);
Cd = app2data(:,3);
AoA = app2data(:,1);
AoAFit = AoA(6:26,1);


Const2 = polyfit(AoA(10:20),Cl(10:20,1),1);
x2 = -10:1:10;
Cl_fit = polyval(Const2,x2);

%% Constants
g = 9.81; %[m/s^2]
AR = 4;
h0 = 7; %[m]

%From milestone 1: Skin Friction Coefficient (Glider/Civil Transport)
Cskin = 0.003;
%Density
rho = 1.14; %[kg/m^3]
%From hand calculations and CAD
Sref = 0.2500; %[m^2]
Swet = 1.4569; %[m^2]
WeightOfFoam = 0.295; %[kgf/m^2]
AoA0 = 0.0775; %2D lift curve slope
AoAL = 0;
e = 0.9;
b = 1; %[m]
chord = 0.25;
S = b*chord;
%fuselage Diameter
Df = 0.1; %[m]
%% Calculations

WeightOfGlider = WeightOfFoam*Swet*g; %[N] Is this right? the foam weight multiplied by the whole area of the glider multiplied by our conversion to newtons?
MassOfGlider = WeightOfGlider/g; %[kg]

%3D CL calculations
CL = AoA0*(AoA-AoAL);
%3D CD
k = 1/(pi*e*AR);
CD = Cd +k*CL.^2;

%Estimating Oswalds
% s = 1-2*(Df/b)^2;
% Q = 1/(.99*s);

% General Milestone 1 calculation
e0 = 1.78*(1-0.045*AR^(0.68))-0.64;
k0 = 1/(pi*e0*AR);

%Parasite Drag
CD0 = Cskin*(Swet/Sref);
% [~,i] = min(abs(CD0-CD));
% CLminD = CL(i);
% CD0 = CD(i) + (k0*CLminD^2);

%Whole Aircraft Drag
CD_Whole = CD0 + (k0*CL.^2);

%L/D
L_D = CL./CD_Whole;
L_DMax = max(L_D);

%Equation for Whole Aircraft Drag Polar
q = polyfit(AoA,CD_Whole,2);


%Angle of Attack
theta = 1/(L_DMax);



%% Max Glide Range
R = h0*L_DMax; %[m]


%% Max Glide Velocity
%from setting CD0=CDi with L/D max
V = sqrt(WeightOfGlider/(sqrt(CD0*(pi*e0*AR))*.5*rho*S)); %[m/s]


%% Max Glide Endurance
E = h0/(V*sin(theta));



%% Results

fprintf('Wing Area: %f [m^2]\n',Sref)
fprintf('Wetted Area: %f [m^2]\n',Swet)
fprintf('Weight: %f [N]\nMass: %f [kg]\n',WeightOfGlider,MassOfGlider)
fprintf('Aspect Ratio [AR]: %f\n',AR)
fprintf('Max Glide Range: %f [m]\n',R)
fprintf('Max Glide Velocity: %f [m/s]\n',V)
fprintf('Max Glide Endurance: %f [s]\n',E)




%% Plot
%Plot Cl
figure(1)
plot(AoA,Cl)
hold on
grid on
plot(AoA,CL)
plot(AoAFit,Cl_fit)
xline(0);
yline(0);
legend('2D Cl', '3D CL','2D Cl fit','x-axis','y-axis','location','southeast')
xlim([-5 12])
title('Coeff. of Lift vs. AoA')
xlabel('Angle of Attack [Degree]')
ylabel('Coeff. of Lift')
hold off


%Plot Cd
figure(2)
plot(AoA,Cd)
hold on
grid on
plot(AoA,CD)
plot(AoA,CD_Whole)
xline(0);
yline(0);
legend('2D Cd', '3D CD','Whole Aircraft CD','x-axis','y-axis','location','southeast')
xlim([-12 12])
title('Coeff. of Drag vs. AoA')
xlabel('Angle of Attack [Degree]')
ylabel('Coeff. of Drag')
hold off
