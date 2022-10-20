%%
clc; clear; close all;

%% Data for just the tarp

[Def1, Load1, ~] = readData("DIC_test1.csv");

L_0 = 12*0.0254; % Initial length [m]
W_0 = 0.75*0.0254; % initial width [m]
t = 0.01*0.0254; % Thickness [m]

Strain1 = Def1/L_0;

A = W_0*t;

Stress1 = Load1/A;



lin_range = 1:174;
Lin_strn1 = Strain1(lin_range);
Lin_strs1 = Stress1(lin_range);

lin_range2 = 174:653;
Lin_strn2 = Strain1(lin_range2);
Lin_strs2 = Stress1(lin_range2);

P11 = polyfit(Lin_strn1, Lin_strs1, 1);
P12 = polyfit(Lin_strn2, Lin_strs2, 1);

LBF1 = polyval(P11,Lin_strn1);
LBF2 = polyval(P12,Lin_strn2);
figure()
scatter(Strain1, Stress1, '.');
hold on
plot(Lin_strn1, LBF1, 'LineWidth', 2);
plot(Lin_strn2, LBF2, 'LineWidth', 2);

%% Finding a polyfit for the entire curve
PFull = polyfit(Strain1(1:end-10),Stress1(1:end-10),9);
FullVal = polyval(PFull,(0:0.001:0.2))';
plot((0:0.001:0.2),FullVal, 'LineWidth', 2)

%% Yield Strength
sig_YS1 = P11(1)*0.02 + P11(2);
plot([0.002 0.002], [0 sig_YS1], 'r--')

sig_max = P11(1)*0.00388 + P11(2);

xlabel("Strain");
ylabel("Stress [Pa]");
title(strcat("Stress v Strain - Young's Modulus of ", num2str(round(P12(1))), " Pa"))

tau_max = max(Stress1);

FOS = sig_YS1/sig_max;

%% Testing with other file

[Def2, Load2, Time2] = readData("Sewn_test1.csv");

Strain2 = Def2/L_0;

A = W_0*t;

Stress2 = Load2/A;

figure;
scatter(Strain2, Stress2, '.');

lin_range = 1:218;
Lin_strn1 = Strain2(lin_range);
Lin_strs1 = Stress2(lin_range);

lin_range2 = 218:653;
Lin_strn2 = Strain2(lin_range2);
Lin_strs2 = Stress2(lin_range2);

P1 = polyfit(Lin_strn1, Lin_strs1, 1);
P2 = polyfit(Lin_strn2, Lin_strs2, 1);

LBF1 = polyval(P1,Lin_strn1);
LBF2 = polyval(P2,Lin_strn2);

hold on;
plot(Lin_strn1, LBF1);
plot(Lin_strn2, LBF2);


%%

function [Def, Load, Time] = readData(flnm)
    data = table2array(readtable(flnm));
    % Crosshead(mm)     Load (N)    Time(s)
    
    Def = data(:,1)/1000; % Deformation, m
    Load = data(:,2); % Load, N
    Time = data(:,3); % Time(s)
end


