%% ASEN 3112: Structures
% Lab 3 - Vibrations
% Author: Nick Bottieri
% Created on: April 11, 2022
% Last updated: April 17, 2022
% Purpose: Plots response on the system as a function of the excitation
% frequency (with a magnification factor applied).

%% Housekeeping
clc
clear all
close all

%% Reading Data file, Manipulating, and Plotting

%%% Read in data:
Data1a = readtable('ASEN3112_Lab3_1');
Data1 = table2array(Data1a);

%%% Fix time data incase the machine time was never reset:
Data1(:,1) = Data1(:,1) - Data1(1,1);

%%% Sampling Frequency + range of frequencies:
sz = size(Data1);
Fs = 1/(Data1(2,1)-Data1(1,1)); % Sampling frequency
f = Fs*(0:(sz(1)/2))/sz(1); % Range of frequencies

%%% FFTs:
tail_fft = fft(Data1(:,3));
wing_fft = fft(Data1(:,4));
nose_fft = fft(Data1(:,5));

%%% Tail fft manipulation:
P2 = abs(tail_fft/sz(1));
P1 = P2(1:sz(1)/2+1);
P1(2:end-1) = 2*P1(2:end-1);

%%% Wing fft manipulation:
P4 = abs(wing_fft/sz(1));
P3 = P4(1:sz(1)/2+1);
P3(2:end-1) = 2*P3(2:end-1);

%%% Nose fft manipulation:
P6 = abs(nose_fft/sz(1));
P5 = P6(1:sz(1)/2+1);
P5(2:end-1) = 2*P5(2:end-1);

%%% Plotting FFT:
figure
hold on
plot(f,P1/max(P1))
plot(f,P3/max(P3))
plot(f,P5/max(P5))
xlim([1 51])
xlabel('Frequency [Hz]','FontSize',16)
ylabel('FFT Amp.','FontSize',16)
title('FFT','FontSize',16)
legend('Tail','Wing','Nose','Location','northwest','FontSize',16)

%%% Mode Shapes
L = 12; % [in]
ne = 2;
nsub = 10;
scale = 1;
ev1 = [12.37*2*pi; 0; 0; 0];
ev2 = [24.50*2*pi; 0; 0; 0];
ev3 = [24.54*2*pi; 0; 0; 0];
ev4 = [42.09*2*pi; 0; 0; 0];
ev5 = [45.67*2*pi; 0; 0; 0];
[x1,v1] = eigenmode(L,ev1,ne,nsub,scale);
[x2,v2] = eigenmode(L,ev2,ne,nsub,scale);
[x3,v3] = eigenmode(L,ev3,ne,nsub,scale);
[x4,v4] = eigenmode(L,ev4,ne,nsub,scale);
[x5,v5] = eigenmode(L,ev5,ne,nsub,scale);

figure
plot(x1,v1)

figure
plot(x2,v2)

figure
plot(x3,v3)

figure
plot(x4,v4)

figure
plot(x5,v5)
