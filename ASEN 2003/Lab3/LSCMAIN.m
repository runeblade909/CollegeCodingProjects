%House Keeping
clc
clear all
close all

%% Intial Variables

r = 7.5; %[cm]
d = 15.5; %[cm]
l = 26; %[cm]
theta0 = 152.5; %[deg]
%% Data Load
%% Column: Time, Wheel Position(deg), Slide Position(mm), Wheel Speed(deg/s), Slide Speed(mm/s), Actual Sample Time (ms)  

[Theta0_5pt,Omega_5pt,Vel_5pt,t_5pt] = LCSDATA('Test1_5pt5v');
[Theta0_6pt,Omega_6pt,Vel_6pt,t_6pt] = LCSDATA('Test1_6pt5v');
[Theta0_7pt,Omega_7pt,Vel_7pt,t_7pt] = LCSDATA('Test1_7pt5v');
[Theta0_8pt,Omega_8pt,Vel_8pt,t_8pt] = LCSDATA('Test1_8pt5v');
[Theta0_9pt,Omega_9pt,Vel_9pt,t_9pt] = LCSDATA('Test1_9pt5v');
[Theta0_10pt,Omega_10pt,Vel_10pt,t_10pt] = LCSDATA('Test1_10pt5v');




%% Velocity Calculations

%if you want to specifiy the number of rotations,


[beta_5,theta_5,velocity_5] = LCSMODEL(r , d, l, Theta0_5pt, Omega_5pt);
[beta_6,theta_6,velocity_6] = LCSMODEL(r , d, l, Theta0_6pt, Omega_6pt);
[beta_7,theta_7,velocity_7] = LCSMODEL(r , d, l, Theta0_7pt, Omega_7pt);
[beta_8,theta_8,velocity_8] = LCSMODEL(r , d, l, Theta0_8pt, Omega_8pt);
[beta_9,theta_9,velocity_9] = LCSMODEL(r , d, l, Theta0_9pt, Omega_9pt);
[beta_10,theta_10,velocity_10] = LCSMODEL(r , d, l, Theta0_10pt, Omega_10pt);

%% Res Calculations

[Res_5, Res_mean5, Res_std5, Res2_5, time5] = Residual(Vel_5pt,velocity_5,t_5pt);
[Res_6, Res_mean6, Res_std6, Res2_6, time6] = Residual(Vel_6pt,velocity_6,t_6pt);
[Res_7, Res_mean7, Res_std7, Res2_7, time7] = Residual(Vel_7pt,velocity_7,t_7pt);
[Res_8, Res_mean8, Res_std8, Res2_8, time8] = Residual(Vel_8pt,velocity_8,t_8pt);
[Res_9, Res_mean9, Res_std9, Res2_9, time9] = Residual(Vel_9pt,velocity_9,t_9pt);
[Res_10, Res_mean10, Res_std10, Res2_10, time10] = Residual(Vel_10pt,velocity_10,t_10pt);

%% Plot Velocity Data

subplot(2,3,1)
plot(Theta0_5pt(1:346),velocity_5(1:346))
hold on
scatter(Theta0_5pt(1:346),Vel_5pt(1:346))
legend('Calculated','Experimental')
title('5v Velocity')
xlabel('Theta [deg]')
ylabel('Velocity [cm/s]')
xlim([0 2600])
ylim([-250 250])
hold off

subplot(2,3,2)
plot(Theta0_6pt(1:216),velocity_6(1:216))
hold on
scatter(Theta0_6pt(1:216),Vel_6pt(1:216))
legend('Calculated','Experimental')
title('6v Velocity')
xlabel('Theta [deg]')
ylabel('Velocity [cm/s]')
xlim([0 2600])
ylim([-250 250])
hold off

subplot(2,3,3)
plot(Theta0_7pt(1:176),velocity_7(1:176))
hold on
scatter(Theta0_7pt(1:176),Vel_7pt(1:176))
legend('Calculated','Experimental')
title('7v Velocity')
xlabel('Theta [deg]')
ylabel('Velocity [cm/s]')
xlim([0 2600])
ylim([-250 250])
hold off

subplot(2,3,4)
plot(Theta0_8pt(1:132),velocity_8(1:132))
hold on
scatter(Theta0_8pt(1:132),Vel_8pt(1:132))
legend('Calculated','Experimental')
title('8v Velocity')
xlabel('Theta [deg]')
ylabel('Velocity [cm/s]')
xlim([0 2600])
ylim([-250 250])
hold off

subplot(2,3,5)
plot(Theta0_9pt(1:109),velocity_9(1:109))
hold on
scatter(Theta0_9pt(1:109),Vel_9pt(1:109))
legend('Calculated','Experimental')
title('9v Velocity')
xlabel('Theta [deg]')
ylabel('Velocity [cm/s]')
xlim([0 2600])
ylim([-250 250])
hold off

subplot(2,3,6)
plot(Theta0_10pt(1:97),velocity_10(1:97))
hold on
scatter(Theta0_10pt(1:97),Vel_10pt(1:97))
legend('Calculated','Experimental')
title('10v Velocity')
xlabel('Theta [deg]')
ylabel('Velocity [cm/s]')
xlim([0 2600])
ylim([-250 250])
hold off

%% Plot Residual Data

figure(2)

subplot(2,3,1)
plot(t_5pt,Res_5)
hold on
title('5v Residual')
xlabel('Time [s]')
ylabel('Residual [cm/s]')
ylim([0 25])
hold off

subplot(2,3,2)
plot(t_6pt,Res_6)
hold on
title('6v Residual')
xlabel('Time [s]')
ylabel('Residual [cm/s]')
ylim([0 25])
hold off

subplot(2,3,3)
plot(t_7pt,Res_7)
hold on
title('7v Residual')
xlabel('Time [s]')
ylabel('Residual [cm/s]')
ylim([0 25])
hold off

subplot(2,3,4)
plot(t_8pt,Res_8)
hold on
title('8v Residual')
xlabel('Time [s]')
ylabel('Residual [cm/s]')
ylim([0 25])
hold off

subplot(2,3,5)
plot(t_9pt,Res_9)
hold on
title('9v Residual')
xlabel('Time [s]')
ylabel('Residual [cm/s]')
ylim([0 25])
hold off

subplot(2,3,6)
plot(t_10pt,Res_10)
hold on
title('10v Residual')
xlabel('Time [s]')
ylabel('Residual [cm/s]')
ylim([0 25])
hold off

figure(3)

Volts = 5:10;
Means = [Res_mean5, Res_mean6, Res_mean7, Res_mean8, Res_mean9, Res_mean10];
STDS  = [Res_std5, Res_std6, Res_std7, Res_std8, Res_std9, Res_std10]
STD_pos = [Res_mean5, Res_mean6, Res_mean7, Res_mean8, Res_mean9, Res_mean10] + STDS;
STD_neg = [Res_mean5, Res_mean6, Res_mean7, Res_mean8, Res_mean9, Res_mean10] - STDS;
plot(Volts,Means)
title('Mean Residual vs Voltage')
xlabel('Voltage [Volts]')
ylabel('Residual [cm/s]')
hold on
plot(Volts,STD_pos,'--r')
plot(Volts,STD_neg,'--r')
legend ('Residuals','Standard Deviation')
hold off

figure(4)

subplot(2,3,1)
plot(time5,Res2_5)
hold on
title('5v Residual')
xlabel('Time [s]')
ylabel('Residual [cm/s]')
ylim([0 25])
hold off

subplot(2,3,2)
plot(time6,Res2_6)
hold on
title('6v Residual')
xlabel('Time [s]')
ylabel('Residual [cm/s]')
ylim([0 25])
hold off

subplot(2,3,3)
plot(time7,Res2_7)
hold on
title('7v Residual')
xlabel('Time [s]')
ylabel('Residual [cm/s]')
ylim([0 25])
hold off

subplot(2,3,4)
plot(time8,Res2_8)
hold on
title('8v Residual')
xlabel('Time [s]')
ylabel('Residual [cm/s]')
ylim([0 25])
hold off

subplot(2,3,5)
plot(time9,Res2_9)
hold on
title('9v Residual')
xlabel('Time [s]')
ylabel('Residual [cm/s]')
ylim([0 25])
hold off

subplot(2,3,6)
plot(time10,Res2_10)
hold on
title('10v Residual')
xlabel('Time [s]')
ylabel('Residual [cm/s]')
ylim([0 25])
hold off