clc; clear all; close all;

%% Vertical Earth Pressure Modeling
% Model found using research paper on vert earth pressure
% Link: https://www.hindawi.com/journals/mpe/2019/8257157/

E = 320; % Modulus of elasticity of soil (MPa)
gamma = 18850; % Volume weight of soil (N/m^3)
D = .25; % Diameter of tunnel (m)
syms h % depth variable to solve for
k0 = .95; % Constant for soil (found in plots in research paper)
theta = 45; % Friction slope angle for soil (max approximate)
k1(h) = (-.015*log(E)+.0133)*log(h/D)+1; % First constant factor
k2(h) = (-.0226*D+1.43)*(h/D)^.1; % Second constant factor
k3(h) = (-.0038*theta+.0524)*log(h/D)-.0902*tand(theta)+1; % 3rd factor
q(h) = k0*k1*k2*k3*gamma*h; % soil pressure wrt depth
%for max pressure
k1max = (-.015*log(E)+.0133)*log(3.11/D)+1;
k2max = (-.0226*D+1.43)*(3.11/D)^.1;
k3max = (-.0038*theta+.0524)*log(3.11/D)-.0902*tand(theta)+1;
qMax = k0*k1max*k2max*k3max*gamma*3.11/1000; % max vertical pressure
disp('Max Pressure(kPa):')
disp(qMax);

figure(1)
hold on
fplot(q/1000,[0 5],'k','HandleVisibility','off')
p = xline(3.11,'r--');
legend([p],'Maximum Predicted Depth')
title('Soil Pressure vs. Depth')
xlabel('Depth (m)')
ylabel('Vertical Earth Pressure (kPa)')
hold off

%% Friction Modeling
L = 4.208;
R = 44;
SoilForce = q*(pi*(D/2)*L);
m = 675;
g = 9.8;
W = m*g;
f_r = .3;
alpha = 21.67;
Angle = acosd((h/R)+sind(90-alpha));
SoilFriction = SoilForce*f_r;
WeightFriction = f_r*W*sind(90-Angle);
TotalFriction = SoilFriction+WeightFriction;

figure(2)
hold on 
fplot(TotalFriction/1000,[0 5],'k','HandleVisibility','off')
p = xline(3.11,'r--');
legend([p],'Maximum Predicted Depth','Location','southeast')
title('Friction Force vs. Depth')
xlabel('Depth (m)')
ylabel('Friction Force (kN)')
hold off

%% Ring Spacing model
syms f % Load force variable
% Stress pressure vs. load force model, slope of .0421 found from data
% points of stress given a load on 1/2" thick .5m diamter ring
Pressure = 0.0421*f; % MPa

% Find load when stress pressure is equal to yield strength of steel (250
% MPa) including FoS
OptimalLoad = double(solve(Pressure>=250))/1.5

Spacing(h) = OptimalLoad/(q*.5*1); % Divide optimal load by Load over area due to pressure
figure(3)
fplot(Spacing,[.1 10],'HandleVisibility','off')
xline(3.11,'r--')
legend('Max Predicted Depth','Location','best')
title('Ring Spacing vs. Depth')
xlabel('Depth (m)')
ylabel('Spacing (m)')


