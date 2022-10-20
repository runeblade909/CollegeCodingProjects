%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE CHALLENGE 9 - Guided Template Script
%
% The purpose of this challenge is to propagate an orbit in a two body
% system for one period, and to plot it's specific energy over time.
%
% To complete the challenge, execute the following steps:
% 1) Set an initial condition vector
% 2) Propagate for exactly period of the orbit
% 3) Calculate the specific energy of the s/c vs. time
% 4) Plot the trajectory, include points for where the trajectory starts,
% ends, and the where the Earth is.
% 5) Plot the change in specific energy vs. time
%
% NOTE: DO NOT change any variable names already present in the code.
% 
% Upload your team's script to Gradescope when complete.
% 
% NAME YOUR FILE AS Challenge9_Sec{section number}_Group{group breakout #}.m 
% ***Section numbers are 1 or 2*** 
% EX File Name: Challenge9_Sec1_Group15.m 
% Group 27
% STUDENT TEAMMATES
% 1) Zak Reichenbach
% 2) Nick Bottieri
% 3) Jack Phillip Davis
% 4) Sydney Ann Walthall
% 5) Elijah Stubbs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Housekeeping
clear variables; close all; clc;


%% Set up
mu =398600.4415; % GM of Earth [km^3/s^2]
r = [10000,0,1000]' ; % initial r vetor [km]
v = [0,7.5574,0]'; % initial v vetor [km/s] %Mass of boddy [kg]
m = 20;          %Mass of smaller body [kg]
a = -mu*(norm(v)^2-2*(mu)/ norm(r))^(-1); % calculating a [km]
T = 2*pi*sqrt(a^3/(mu)); % calculating T [s]
y0 = [r;v]; % initial condition vector
y = [r;v];
tspan = [0 T]; % time domain [s] One orbital period

% Propagate w/ ode45
[t,y] = ode45(@(t,y) g_fun_orbit(r,y,mu,t), tspan, y0);

MagV = v(2);           %Starting velocity
Radius = mu/MagV^2;    %Radius between objects
UnitVec = r/norm(r);   %The unit vector between masses
MagEarthRadius = norm(r) - Radius; % Finds the radius from the origin to the earth 
EarthRadius = MagEarthRadius*UnitVec;      %Puts the radius from earth on the unit vector

velo = zeros(121,1);   %Preallocating vectors for energy 
radi = zeros(121,1);
energy = zeros(121,1);
for i =1:121
%Calculates combined vector positions and velocities for Specific Energy
velo(i) = norm(y(i,4:6));
radi(i) = norm(y(i,1:3));
energy(i) = ((velo(i)^2)/2)-(mu/radi(i));
end
%% Plotting
figure(1)
plot3(y(1,1),y(1,2),y(1,3),'g*');           % plot starting point
label1 = {'Start Point'};
text(y(1,1),y(1,2),y(1,3),label1);
hold on; grid minor;
xlabel('x [km]'); 
ylabel('y [km]');
zlabel('z [km]');

title('Orbital Path Of A Satellite');

plot3(y(:,1),y(:,2),y(:,3));                % plot trajectory
plot3(y(121,1),y(121,2),y(121,3),'r*');     % plot ending point
label2 = {'End Point'};
text(y(121,1),y(121,2),y(121,3),label2);

plot3(EarthRadius(1),EarthRadius(2),EarthRadius(3),'bo'); % plot earth
label3 = {'Earth'};
text(EarthRadius(1),EarthRadius(2),EarthRadius(3),label3);

figure(2)
plot(t,energy); %plot specific energy vs. time
grid minor;
xlabel('Time[Sec]');
ylabel('Specific Energy [J/kg]');
title('Specific Energy Vs. Time');


function dr_dt = g_fun_orbit(r,y,mu,t)

%y1 = Rx
%y2 = Ry
%y3 = Rz
%y4 = Vx
%y5 = Vy
%y6 = Vz

y1 = y(1);%position
y2 = y(2);%position
y3 = y(3);%position
y4 = y(4);%Velocity
y5 = y(5);%Velocity
y6 = y(6);%Velocity

y1_dot = y4;
y2_dot = y5;
y3_dot = y6;
y4_dot = (-mu/norm(r)^3)*y1;
y5_dot = (-mu/norm(r)^3)*y2;
y6_dot = (-mu/norm(r)^3)*y3;


dr_dt = [y1_dot;y2_dot;y3_dot;y4_dot;y5_dot;y6_dot];


end
