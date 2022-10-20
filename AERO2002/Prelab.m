clear all
clc
%% Gathering and Converting Data
data = readtable('Balloon Prelab Data 2020.xlsx', 'FileType', 'spreadsheet');
%Importing Data
data = data(1:end-1, :);
%Mass Balloon
ballmasskg = mean(data.BalloonMass_g_(1:3)./1000);
sigballmasskg = data.BalloonMass_g_(4:end)./1000;
%Mass Sting
stringmasskg = mean(data.StringMass_g_(1:3)./1000);
sigstringmasskg = data.StringMass_g_(4:end)./1000;
%Mass Added Mass
addmasskg = mean(data.AdditionalMass_g_(1:3)./1000);
sigaddmasskg = data.AdditionalMass_g_(4:end)./1000;
%Payload for Neutral Buoyancy
paymassg = mean(str2double(data.PayloadMass__forNeutralBuoyancy__g_(1:3)));
sigpaymassg = str2double(data.PayloadMass__forNeutralBuoyancy__g_(4:end));

paymasskg = mean(paymassg./1000);
sigpaymasskg = sigpaymassg./1000;
%Tank Height
waterHeightm = mean(str2double(data.WaterHeightDisplacement_mm_(1:3))./1000); %[m]
sigwaterHeightm = str2double(data.WaterHeightDisplacement_mm_(4:end))./1000; %[m]
%Conversion
in2m = 0.0254;
%% Pressure
P = 840*100; %[Pa]
sigP = .05; %[Pa]
%% Temperature
T = (71-32)*(5/9)+273.15; %[F]
sigT = .5; %[F]
%% Diameter and Thickness
BallD = 13.5; %[in]
sigBallD = 1/16; %[in]
BallDm = BallD * in2m;
BallT = 8.5; %[in]
BallTm = BallT * in2m;
% v = pi * r^2 * thickness
MesuredVol = pi*(BallDm/2)^2 * BallTm;
%% Converting to Meters
tankLin = 19.5; %[in]
sigtankLin = 1/16; %[in]
tankWin = 9.75; %[in]
sigtankWin = 1/16; %[in]
tankL = tankLin*in2m; %[M]
tankW = tankWin*in2m; % [M]
sigtankL = sigtankLin*in2m; %[M]
sigtankW = sigtankWin*in2m; % [M]

%% Math Section
%AirDensity = 1.14; % kg/m^3
R = 8.3145; %J/M/K
molarmassAirKg = 28.97/1000;
molarmassHe = 4.002602; %[kg]
BallVol = waterHeightm.*tankL.*tankW; %[m^3]
DensityOfHelium = .164; %[kg/m^3]
%massOfAirDisplaced = 1.292*BallVol; %estimated air density at 
n = P*BallVol/(R*T);
%massOfAirDisplaced = n*molarmassAirKg;
totalMass = mean(ballmasskg+stringmasskg+paymasskg);
AirDensity = P/(R/molarmassAirKg)/T;

%Force Balance roe*g*V = Mtot * g % YOU DONT NEED g, DONT WORRY ABOUT IT
%Mtot = totalMass + MassHe
Fb = AirDensity*BallVol %* 9.8;
FGasWeight = Fb - totalMass %*(9.8);

%expect 1 or 2 grams
