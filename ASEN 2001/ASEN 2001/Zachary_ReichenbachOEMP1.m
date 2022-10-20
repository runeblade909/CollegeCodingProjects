clc
clear all
%% Load Data
load('OEMP1_data.mat')


%% Graph
 MachData = table2array(flight_envelope(1:end,1));
 AltitudeData = (table2array(flight_envelope(1:end,2))*3.28084)/1000;%in feet
 AirDensityMet = table2array(flight_envelope(1:end,4)); %kg/m^3
 AirDensityImp = AirDensityMet.*(2.204621)*(1/3.28084)^3; %lb/ft^3
%% Force F = A x Q x Cd
%A = 1/2 x b x h
 Area = 1/2*(.4)*(1); % ft^2
 SpeedOfSound = table2array(flight_envelope(:,3))*3.28084;
 FPS = SpeedOfSound;
 Cd = 1.2;
 
%Q = 1/2 x Roe x v^2;

Q =  (1/2) .* AirDensityImp .* FPS.^2; %lb/ft^3 * (ft/s)^2

F = Area * Q * Cd; %(lb*ft)/s^2 or SLUGS 

%% Graph
figure(1)
plot(MachData,AltitudeData);
ylabel('Height (kft) aircraft is flying at')
xlabel('Speed of aircraft (Mach)')
title('Sustained Flight Envelope for F-105')

figure(2)
plot(AltitudeData,F);
ylabel('Force (slugs)')
xlabel('Height (kft) aircraft is flying at')
title('Worst Case Scenario Point Load vs Altitude')
hold on
%% Maximum point load
[MaxLoad,i] = max(F);
fprintf('The maximum load that the ram jet blade faces is %.2f slugs.\nThis occured at an altitude of %.2f feet,\nand a speed of %.2f feet per second.',MaxLoad,AltitudeData(i)*1000,FPS(i));
xline(AltitudeData(i),'-','Maximum');
hold off




