clear;clc;close all;

%Code by: Zhixing Yao
%Edited and Analysed by: Zak Reichenbach

%% Experimental Data
path = "Rigid Arm Data/";
files = dir(path+"*.txt");
filename = files(randi([1,length(files)])).name;
data = load(path+filename);
time = (data(:,1)-data(1,1))./1000;
theta = data(:,2);

%% Constants
Kg = 33.3; % gear ratio
Km = 0.0401; % motor constant
Rm = 19.2; % armature resistance [ohm]
J = 0.0005+0.0015; % total inertia [kg*m^2]

fpara = strsplit(filename,"_");
Kpt = str2double(fpara{1});
Kdt = str2double(fpara{2});

%% Transfer Function
num = Kpt*Kg*Km/(J*Rm);
den = [1,Kg^2*Km^2/(J*Rm)+Kdt*Kg*Km/(J*Rm),Kpt*Kg*Km/(J*Rm)];
sysTF = tf(num,den);

s = stepinfo(sysTF,'SettlingTimeThreshold',.05);

T = time;
U = data(:,6);
[thetaSim,timeSim,~] = lsim(sysTF,U,T);
hold on;
plot(timeSim,thetaSim,"LineWidth",2);
plot(time,theta,"LineWidth",2);
title("Rigid Arm Kpt: "+Kpt+"  Kdt: "+Kdt);
legend("simulation","experiment");
ylabel("theta [rad]");
xlabel("time [s]");