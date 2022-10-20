%% Header
% Zak Reichenbach
% Elena Bauer
% Matt Gozum
% Gavin O'Connell
%ASEN 3128
%Lab_1.m
%Created: 8/24/21

%% Housekeeping
clear
clc
close all
%% Problem 1

tspan = [0 5];  %Timespan [seconds]
y0 = [8; 5; 2]; % [Initial Conditions]
[t,f] = ode45(@(t,f) ODEfun(t,y0), tspan, y0);

%Plots
figure('NumberTitle', 'off', 'Name', '1.a Basic ODE')
title('1.a')
subplot(3,1,1)
plot(t,f(:,1),'r', 'LineWidth', 2)
title('xdot')
subplot(3,1,2)
plot(t,f(:,2),'g', 'LineWidth', 2)
title('ydot')
subplot(3,1,3)
plot(t,f(:,3),'b', 'LineWidth', 2)
title('zdot')

%% Problem 2

%Constants
m = 0.03;               %kg
dia = 0.03;             %m
A = pi*(dia^2)/4;       %m^2
rho = 1.225;            %kg/m^3 at S.L.
Cd = 0.6;               %Coeff. of Drag
g = 9.81;               %m/s^2
x0 = [0,0,0,0,20,-20]'; %Initial Conditions
wind = [0,0,0];         %Preallocated 3D Wind fector
tspan = [0 200];        %Time Span [seconds]
opt = odeset('Events', @groundStrike);
[t1,j] = ode45(@(t,x) objectEOM(t,x,rho,Cd,A,m,g,wind), tspan, x0, opt);

 dest = [j(end,1), j(end,2)]; %Landing Location [x,y]
 distance = norm(dest);       %Horizontal Distance [m]


%Plots
figure()
hold on
plot3(j(:,1), j(:,2), -j(:,3))
plot3(j(1,1),j(1,2),-j(1,3),'go')
zlabel('Z axis (flipped)')
ylabel('y axis')
xlabel('x axis')
grid on
hold off

%% Problem 2b
N = 100;                    %Number of test winds
randwind = rand(3,N)*100;   %Wind values
randwind(2:3,:) = 0;        %Horizontal Wind only    
distance2 = zeros(1,N);     %Preallocated wind array

trials = struct()

for i =1 : N
    
    wind = randwind(:,i);   
    
    windmag(i) = norm(wind);
    
    
    
    %j2 = [x y z Vx Vy Vz]
    
    [t1,j2] = ode45(@(t,x) objectEOM(t,x,rho,Cd,A,m,g,wind), tspan, x0, opt);
    
   
    dest2 = [j2(end,1), j2(end,2)];
    distance2(i) = norm(dest2);
        
end

figure()
    hold on
    grid on
    plot3(j2(:,1), j2(:,2), -j2(:,3))
    grid on
hold off

NoWindVec = [distance;0;0;0;distance];

Comparison = [randwind; windmag ; distance2];

Comparison = [NoWindVec Comparison];

figure()

scatter(Comparison(5,:), Comparison(4,:));
xlabel("Wind Speed Magnitude [m/s]")
ylabel("Distance [m]")

%% Function Junction
function xdot = objectEOM(t,x,rho,Cd,A,m,g,wind)
xdot(1) = x(4);
xdot(2) = x(5);
xdot(3) = x(6);

airrel = x(4:6)- wind;
airspeed = norm(airrel);
DragMag = Cd*(rho*(airspeed)^2)*0.5*A;
unitvec = airrel/airspeed;
DragF = -DragMag * unitvec;


xdot(4) = DragF(1)/m;
xdot(5) = DragF(2)/m;
xdot(6) = DragF(3)/m + g;
xdot = xdot';
end

function der = ODEfun(t,y0)
x = y0(1);
y = y0(2);
z = y0(3);
dxdt = x + 2*y + z;
dydt = x - 5*z;
dzdt = x*y - y.^2 + 3*z.^3;
der = [dxdt; dydt; dzdt];
end

% ground strike func
function [value, isterminal, direction] = groundStrike(t,x)
    value = (x(3) <= 0);
    isterminal = 1;
    direction = 0;
end