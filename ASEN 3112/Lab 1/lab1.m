% Ashwin Raju
% ASEN 3112
% lab1.m
% Created: 2/04/2022

%% Housekeeping
clear
clc
close all

%% Constants
De = .75; % Exterior Diameter [in]
thic = 1/16; % Uniform Wall Thickness [in]
L = 1; % Extensometer Gage Length [in]
distGrips = 10; % Distance between Grips [in]
G = 3.75 * 10^6; % Shear Modulus [psi]
Re = De/2; % Exterior Tube Radius [in]

Re = (Re +(Re-thic))/2;


%% Closed Wall Data

% Import and Order Data
closedWallData = load('400inlbf_03.txt');
t_Closed = closedWallData(:,1); %time [s]
twistAngle_Closed = deg2rad( closedWallData(:,2) ); % Twist Angle [rad]
gammaExt_Closed = deg2rad( closedWallData(:,3) ); % Extensometer Strain [rad]
torque_Closed = closedWallData(:,4); % Torque [in*lbs]

%Turn data
deltaRadC = max(closedWallData(:,2)) - min(closedWallData(:,2));
ShearMaxC = max(abs(closedWallData(:,3)));


% Zero the angle data
twistAngle_Closed = twistAngle_Closed - twistAngle_Closed(1);
gammaExt_Closed = gammaExt_Closed - gammaExt_Closed(1);

% Find strain from twist angle data
gammaGrip_Closed = twistAngle_Closed*(Re/distGrips); 


% Line Fitting Extensometer Strain
line1 = fitlm(gammaExt_Closed, torque_Closed, 'poly1');
coeffs1 = table2array(line1.Coefficients);
lineExt_Closed = [coeffs1(2,1) coeffs1(1,1)];
yExt_Closed = polyval(lineExt_Closed, gammaExt_Closed);
GJ_Ext_Closed = lineExt_Closed(1)*Re;
sigmaGJ1 = coeffs1(2,2)*Re;

% Line Fitting Encoder Strain
line2 = fitlm(gammaGrip_Closed, torque_Closed, 'poly1');
coeffs2 = table2array(line2.Coefficients);
lineGrip_Closed = [coeffs2(2,1) coeffs2(1,1)];
yGrip_Closed = polyval(lineGrip_Closed, gammaGrip_Closed);
GJ_Grip_Closed = lineGrip_Closed(1)*Re;
sigmaGJ2 = coeffs2(2,2)*Re;

%Zaks Stuff

TmaxC = max(abs(torque_Closed));

lineGrip_Closed(1);

SlopeC = TmaxC/ShearMaxC;

GJC = SlopeC*Re;


GGJJC = torque_Closed./twistAngle_Closed;
GGJJC = mean(GGJJC(2:end))



%% Plots
figure(1);
subplot(2,1,1); 
plot(gammaExt_Closed, torque_Closed,'LineWidth',1.5);
hold on;
title('Extensiometer \gamma vs. Torque');
ylabel('Torque [in*lbs]');
xlabel('\gamma [rad]');
plot(gammaExt_Closed, yExt_Closed,'LineWidth',1.5)
GJstr1 = sprintf('GJ = %.1f',GJ_Ext_Closed);
legend('Torque-Strain Curve', GJstr1,'Location','Best')

subplot(2,1,2);
plot(gammaGrip_Closed, torque_Closed,'LineWidth',1.5);
hold on
title('Calculated Encoder \gamma vs. Torque');
ylabel('Torque [in*lbs]');
xlabel('\gamma [rad]');
plot(gammaGrip_Closed, yGrip_Closed,'LineWidth',1.5);
GJstr2 = sprintf('GJ = %.1f',GJ_Grip_Closed);
legend('Torque-Strain Curve', GJstr2,'Location','Best')

sgtitle('Closed Thin Wall Specimen') 
%% Open Wall Data

% Import and Order Data
openWallData = load('20inlbf_03.txt');
t_Open = openWallData(:,1); %time [s]
twistAngle_Open = deg2rad( openWallData(:,2) ); % Twist Angle [rad]
gammaExt_Open = deg2rad( openWallData(:,3) ); % Extensometer Strain [rad]
torque_Open = openWallData(:,4); % Torque [in*lbs]


%Turn data
deltaRadO = max(openWallData(:,2)) - min(openWallData(:,2));
ShearMaxO = max(abs(openWallData(:,3)));


% Zero the angle data
twistAngle_Open = twistAngle_Open - twistAngle_Open(1);
    gammaExt_Open = gammaExt_Open - gammaExt_Open(1);

% Find strain from twist angle data
%This is not the recovery region, the recovery region. the line below with
%steaper slope, uses thic
    gammaGrip_Open = twistAngle_Open*(Re/distGrips);

% Line Fitting Extensometer Strain
line3 = fitlm(gammaExt_Open, torque_Open, 'poly1');
coeffs3 = table2array(line3.Coefficients);
lineExt_Open = [coeffs3(2,1) coeffs3(1,1)];
yExt_Open = polyval(lineExt_Open, gammaExt_Open);
GJ_Ext_Open = lineExt_Open(1)*thic;
sigmaGJ3 = coeffs3(2,2)*thic;

% Line Fitting Encoder Strain
line4 = fitlm(gammaGrip_Open, torque_Open, 'poly1');
coeffs4 = table2array(line4.Coefficients);
lineGrip_Open = [coeffs4(2,1) coeffs4(1,1)];
yGrip_Open = polyval(lineGrip_Open, gammaGrip_Open);
GJ_Grip_Open = lineGrip_Open(1)*thic;
sigmaGJ4 = coeffs4(2,2)*thic;

%Zaks Stuff

TmaxO = max(abs(torque_Open));

lineExt_Open(1);

SlopeO = TmaxO/ShearMaxO;

GJO = SlopeO*Re;

GGJJO = torque_Open./twistAngle_Open;
% GGJJO = mean(GGJJO(2:end))


%% Plots
figure(2);

subplot(2,1,1);
plot(gammaExt_Open,torque_Open);
hold on
title('Extensiometer \gamma vs. Torque');
ylabel('Torque [in*lbs]');
xlabel('\gamma [rad]');
plot(gammaExt_Open, yExt_Open,'LineWidth',1.5)
GJstr3 = sprintf('GJ = %.1f',GJ_Ext_Open);
legend('Torque-Strain Curve', GJstr3,'Location','Best')

subplot(2,1,2); 
plot(gammaGrip_Open,torque_Open);
hold on
title('Calculated Encoder \gamma vs. Torque');
ylabel('Torque [in*lbs]');
xlabel('\gamma [rad]');
plot(gammaGrip_Open, yGrip_Open,'LineWidth',1.5);
GJstr4 = sprintf('GJ = %.1f',GJ_Grip_Open);
legend('Torque-Strain Curve', GJstr4,'Location','Best')



sgtitle('Open Thin Wall Specimen') 













