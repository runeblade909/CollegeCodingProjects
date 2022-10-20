%% Header
%   Author:
%       Cody Wheeler, but mostly Greyson Dudley (The Propulsion God)
%   SensorReading.m
%   Created: 09/08/21
%   Purpose: 

%%% IMPORTANT: YOU NEED "remove_nan.m" TO RUN THIS CODE AS WELL AS TEST DATA A AND B!!!
%% Housekeeping
close all;
clc; clear;

%% Read File

data = load("test_data_A.mat");
data2 = load("test_data_B.mat");

data = data.data_A(:,:);
data2 = data2.data_A(:,:);
for i = 2:length(data)-1
    data(i,5) = (i - 1) +  (2 * 0.0254);
    data2(i,5) = (i - 1) +  (2 * 0.0254);
end

data = remove_nan(data);
data2 = remove_nan(data2);

lat_distance = data(1:102,5); % distance pulled (0 at tunnel end), m
vert_distance1 = data(1:102,1)/1000; % distance from furthest sensor to top of tunnel, m
vert_distance2 = data(1:102,2)/1000; % distance from second sensor to top of tunnel, m
vert_distance3 = data(1:102,3)/1000; % distance from third sensor to top of tunnel, m
vert_distance4 = data(1:102,4)/1000; % distance from fourth sensor to top of tunnel, m


Unloaded_lat = data2(1:102,5); % unloadded distance pulled, m
Unloaded_vert1 = data2(1:102,1)/1000; % unloaded distance for furthest sensor, m
Unloaded_vert2 = data2(1:102,2)/1000; % unloaded distance for furthest sensor, m
Unloaded_vert3 = data2(1:102,3)/1000; % unloaded distance for furthest sensor, m
Unloaded_vert4 = data2(1:102,4)/1000; % unloaded distance for furthest sensor, m

%% Constants
        %%Find distances apart for sensors
Sensor_dist1 = 0/1000;
Sensor_dist2 = 50/1000;
Sensor_dist3 = 100/1000;
Sensor_dist4 = 150/1000;

ring_diameter = .5173; % center diameter, m 
circumference = ring_diameter*pi; % center circumference, m
% E_rod = ; % Youngs modulus of steel rod, Pa
E_tarp = 1.07*10^9; % Youngs Modulus of tarp, Pa
A = 2*pi*(ring_diameter/2)^2; % Surface area, m^2

Max_length = max(lat_distance); % tunnel length, m
Max_depth = 0.8; % max depth of sensor, m
Angle = atan(Max_depth/Max_length); % Tunnel angle, deg

Depth_Sensor1 = 0;
Depth_Sensor2 = (Max_length/(Max_length - Sensor_dist2)) * Max_depth - Max_depth;
Depth_Sensor3 = (Max_length/(Max_length - Sensor_dist3)) * Max_depth - Max_depth;
Depth_Sensor4 = (Max_length/(Max_length - Sensor_dist4)) * Max_depth - Max_depth;

%% Initial data plot and curve fit
p1 = polyfit(lat_distance + Sensor_dist1,vert_distance1,3);
x1 = linspace(Sensor_dist1, Max_length, length(lat_distance));
y1 = polyval(p1,x1);

p2 = polyfit(lat_distance + Sensor_dist2,vert_distance2,3);
x2 = linspace(Sensor_dist2, Max_length, length(lat_distance));
y2 = polyval(p2,x2);

p3 = polyfit(lat_distance + Sensor_dist3,vert_distance3,3);
x3 = linspace(Sensor_dist3, Max_length, length(lat_distance));
y3 = polyval(p3,x3);

p4 = polyfit(lat_distance + Sensor_dist4,vert_distance4,3);
x4 = linspace(Sensor_dist4, Max_length, length(lat_distance));
y4 = polyval(p4,x4);

figure
hold on 
plot(lat_distance + Sensor_dist2,vert_distance2);
plot(x2,y2);
plot(lat_distance + Sensor_dist3,vert_distance3);
plot(x3,y3);
plot(lat_distance + Sensor_dist4,vert_distance4);
plot(x4,y4);
title('Sensor Distance Measurement vs. Distance Along Tunnel')
xlabel('Distance Along Tunnel (m)')
ylabel('Sensor Measurement (m)')
legend('Sensor 1/2','Curve Fit 1/2','Sensor 3','Curve Fit 3','Sensor 4','Curve Fit 4','location','best')
hold off




%% Deformation, stress, and strain analysis
Deformation1 = Unloaded_vert1 - vert_distance1; % Deformation between loaded and unloaded measurement from sensor 1
Deformation2 = Unloaded_vert2 - vert_distance2; % Deformation between loaded and unloaded measurement from sensor 2
Deformation3 = Unloaded_vert3 - vert_distance3; % Deformation between loaded and unloaded measurement from sensor 3
Deformation4 = Unloaded_vert4 - vert_distance4; % Deformation between loaded and unloaded measurement from sensor 4

Strain1 = Deformation1./Unloaded_vert1;
Strain2 = Deformation2./Unloaded_vert2;
Strain3 = Deformation3./Unloaded_vert3;
Strain4 = Deformation4./Unloaded_vert4;

%% SHOULD USE HOOP STRESS, BUT WE NEED PRESSURE
% Hoop stress calculations
d_rod = 0.5*2.54/100; % m
E_rod = 215e9; %GPa
I_rod = pi*d_rod^4/32; % m^4
R_ring = ring_diameter/2;

Stress1 = Deformation1 * E_rod * I_rod / (0.149 * R_ring^3);
Stress2 = Deformation2 * E_rod * I_rod / (0.149 * R_ring^3);
Stress3 = Deformation3 * E_rod * I_rod / (0.149 * R_ring^3);
Stress4 = Deformation4 * E_rod * I_rod / (0.149 * R_ring^3);

% Stress1 = E_tarp*Strain1; % Stress from first sensor data, Pa
% Stress2 = E_tarp*Strain2; % Stress from second sensor data, Pa
% Stress3 = E_tarp*Strain3; % Stress from third sensor data, Pa
% Stress4 = E_tarp*Strain4; % Stress from fourth sensor data, Pa

% Load1 = Stress1*A;
% Load2 = Stress2*A;
% Load3 = Stress3*A;
% Load4 = Stress4*A;

% What does this do?
Load1 = (Stress1*.254/1000)/(ring_diameter/2);
Load2 = (Stress2*.254/1000)/(ring_diameter/2);
Load3 = (Stress3*.254/1000)/(ring_diameter/2);
Load4 = (Stress4*.254/1000)/(ring_diameter/2);

Depth_unloaded = ((Max_length - lat_distance)./Max_length)' .* Max_depth;
Depth_loaded = ((Max_length - Unloaded_lat)./Max_length)' .* Max_depth;

p5 = polyfit(Depth_loaded - Depth_Sensor1,Load1,3);
x5 = linspace(Depth_Sensor1, Max_depth, length(lat_distance));
y5 = polyval(p5,x5);

p6 = polyfit(Depth_loaded - Depth_Sensor2,Load2,3);
x6 = linspace(Depth_Sensor2, Max_depth, length(lat_distance));
y6 = polyval(p6,x6);

p7 = polyfit(Depth_loaded - Depth_Sensor3,Load3,3);
x7 = linspace(Depth_Sensor3, Max_depth, length(lat_distance));
y7 = polyval(p7,x7);

p8 = polyfit(Depth_loaded - Depth_Sensor4,Load4,3);
x8 = linspace(Depth_Sensor4, Max_depth, length(lat_distance));
y8 = polyval(p8,x8);

figure
hold on 
plot(Depth_loaded - Depth_Sensor2,Load2);
plot(x6,y6);
plot(Depth_loaded - Depth_Sensor3,Load3);
plot(x7,y7);
plot(Depth_loaded - Depth_Sensor4,Load4);
plot(x8,y8);
title('Measured Load vs. Depth')
xlabel('Depth (m)')
ylabel('Load (Pa)')
legend('Load Measurement 1/2','Curve Fit 1/2','Load Measurement 3','Curve Fit 3','Load Measurement 4','Curve Fit 4','location','best')
hold off

%% Stress/Pressure Model
p9 = polyfit(Depth_loaded - Depth_Sensor1,Stress1,3);
x9 = linspace(Depth_Sensor1, Max_depth, length(lat_distance));
y9 = polyval(p9,x9);
StressEq_1 = p9(3) + p9(2)*x9.^2 + p9(3)*x9.^3;

p10 = polyfit(Depth_loaded - Depth_Sensor2,Stress2,3);
x10 = linspace(Depth_Sensor2, Max_depth, length(lat_distance));
y10 = polyval(p10,x10);
StressEq_2 = p10(3) + p10(2)*x10.^2 + p10(3)*x10.^3;

p11 = polyfit(Depth_loaded - Depth_Sensor3,Stress3,3);
x11 = linspace(Depth_Sensor3, Max_depth, length(lat_distance));
y11 = polyval(p11,x11);
StressEq_3 = p11(3) + p11(2)*x11.^2 + p11(3)*x11.^3;

p12 = polyfit(Depth_loaded - Depth_Sensor4,Stress4,3);
x12 = linspace(Depth_Sensor4, Max_depth, length(lat_distance));
y12 = polyval(p12,x12);
StressEq_4 = p12(3) + p12(2)*x12.^2 + p12(3)*x12.^3;

figure
hold on 
plot(Depth_loaded - Depth_Sensor2,Stress2);
plot(x10,y10);
plot(Depth_loaded - Depth_Sensor3,Stress3);
plot(x11,y11);
plot(Depth_loaded - Depth_Sensor4,Stress4);
plot(x12,y12);
title('Measured Stress vs. Depth')
xlabel('Depth (m)')
ylabel('Stress (Pa)')
legend('Stress Measurement 1/2','Curve Fit 1/2','Stress Measurement 3','Curve Fit 3','Stress Measurement 4','Curve Fit 4','location','best')
hold off

%% Deformation Model
p13 = polyfit(Depth_loaded - Depth_Sensor1,Deformation1,3);
x13 = linspace(Depth_Sensor1, Max_depth, length(lat_distance));
y13 = polyval(p13,x13);

p14 = polyfit(Depth_loaded - Depth_Sensor2,Deformation2,3);
x14 = linspace(Depth_Sensor2, Max_depth, length(lat_distance));
y14 = polyval(p14,x14);

p15 = polyfit(Depth_loaded - Depth_Sensor3,Deformation3,3);
x15 = linspace(Depth_Sensor3, Max_depth, length(lat_distance));
y15 = polyval(p15,x15);

p16 = polyfit(Depth_loaded - Depth_Sensor4,Deformation4,3);
x16 = linspace(Depth_Sensor4, Max_depth, length(lat_distance));
y16 = polyval(p16,x16)

figure
hold on 
plot(Depth_loaded - Depth_Sensor2,Deformation2);
plot(x14,y14);
plot(Depth_loaded - Depth_Sensor3,Deformation3);
plot(x15,y15);
plot(Depth_loaded - Depth_Sensor4,Deformation4);
plot(x16,y16);
title('Measured Deformation vs. Depth')
xlabel('Depth (m)')
ylabel('Deformation (m)')
legend('Deformation Measurement 1/2','Curve Fit 1/2','Deformation Measurement 3','Curve Fit 3','Deformation Measurement 4','Curve Fit 4','location','best')
hold off

%% Load vs Deformation

figure;
hold on;
scatter(Stress2, Deformation2)
scatter(Stress3, Deformation3)
scatter(Stress4, Deformation4)

legend("Load 1/2", "Load 3", "Load 4");
title("Stress applied vs deformation");
xlabel("Stress (Pa)");
ylabel("Deformation (m)");

