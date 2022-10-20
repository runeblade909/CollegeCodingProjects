
%% ASEN-2003 Lab4: Balanced and Unbalanced Wheel
% Purpose: Load in apparatus data, call model functions, and plot results
% Authors: Jacob Starkel, Krystal Horton, Cate Billings, and Zak Reichenbach
% Date Completed: 3/12/2021

% House Keeping
clc
clear
close all

%% Column: Time, Wheel Position(deg), Slide Position(mm), Wheel Speed(deg/s), Slide Speed(mm/s), Actual Sample Time (ms)
%Pulling the data we need for model calculations

[Theta_b1,Omega0_b1,t_b1] = UBDATA('balanced_1');
[Theta_b2,Omega0_b2,t_b2] = UBDATA('balanced_2');
[Theta_u1,Omega0_u1,t_u1] = UBDATA('unbalanced_1');
[Theta_u2,Omega0_u2,t_u2] = UBDATA('unbalanced_2');

% Sorted by eye for a range of roughly 0.5<theta<15 radians 

Theta0_b1 = Theta_b1(70:106);
Theta0_b2 = Theta_b2(69:106);
Theta0_u1 = Theta_u1(9:37);
Theta0_u2 = Theta_u2(9:36);

Omega_b1 = Omega0_b1(70:106);
Omega_b2 = Omega0_b2(69:106);
Omega_u1 = Omega0_u1(9:37);
Omega_u2 = Omega0_u2(9:36);
%% Balanced Models
%% Model 1 Function Calls

w = model_1(Theta0_b1);

wb = model_1(Theta0_b2);

%% Model 2 Function Calls

w_2 = model_2(Theta0_b1);

w_2b = model_2(Theta0_b2);


%% Model 1 & 2 Report Plots

figure(1)
subplot(2,1,1)
hold on
plot(Theta0_b1,w,'g');
plot(Theta0_b1,w_2(:,1),'kd-','MarkerSize',3);
plot(Theta0_b1,w_2(:,2),'ko-','MarkerSize',3);
plot(Theta0_b1,w_2(:,3),'ks-','MarkerSize',3);
plot(Theta0_b1,w_2(:,4),'kx-','MarkerSize',3);
plot(Theta0_b1,w_2(:,5),'b*-','MarkerSize',3);
scatter(Theta0_b1,Omega_b1,22,'rx')
legend('Model 1','T = 1 [N*m]','T = 1.2 [N*m]','T = 1.4 [N*m]','T = 1.45 [N*m]','Model 2 : T = 1.5 [N*m] ','Experimental Data')
xlabel('Angular Position [rad]')
ylabel('Angular Velocity [rad/s]')
title('Angular Velocity vs Angular Position')
hold off

subplot(2,1,2)
hold on
plot(Theta0_b1,w_2(:,5),'b','MarkerSize',3);
scatter(Theta0_b1,Omega_b1,21,'rx')
scatter(Theta0_b2,Omega_b2,21,'gx')
legend('Model 2','Experimental Data 1','Experimental Data 2')
xlabel('Angular Position [rad]')
ylabel('Angular Velocity [rad/s]')
title('Angular Velocity vs Angular Position')
hold off

figure(5)
hold on
plot(Theta0_b1,w,'g');
plot(Theta0_b1,w_2(:,1),'kd-','MarkerSize',3);
plot(Theta0_b1,w_2(:,2),'ko-','MarkerSize',3);
plot(Theta0_b1,w_2(:,3),'ks-','MarkerSize',3);
plot(Theta0_b1,w_2(:,4),'kx-','MarkerSize',3);
plot(Theta0_b1,w_2(:,5),'b*-','MarkerSize',3);
scatter(Theta0_b1,Omega_b1,22,'rx')
legend('Model 1','T = 1 [N*m]','T = 1.2 [N*m]','T = 1.4 [N*m]','T = 1.45 [N*m]','Model 2 : T = 1.5 [N*m] ','Experimental Data')
xlabel('Angular Position [rad]')
ylabel('Angular Velocity [rad/s]')
title('Angular Velocity vs Angular Position')
hold off
%% Unbalanced Models
%% Model 3 Function Calls
w_3 = model_3(Theta0_u2);

%% Model 4 Function Calls
w_4 = model_4(Theta0_u2);

%% Model 3 and 4 Report Plots
figure(2)
hold on
plot(Theta0_u2,w_3,'k');
plot(Theta0_u2,w_4,'b');
scatter(Theta0_u1,Omega_u1,20,'rx');
scatter(Theta0_u2,Omega_u2,20,'gx');
legend('Model 3','Model 4','Unbalanced Experimental 1','Unbalanced Experimental 2')
xlabel('Angular Position [rad]')
ylabel('Angular Velocity [rad/s]')
title('Angular Velocity vs Angular Position')

figure(6)
hold on
plot(Theta0_b1,w,'k');
plot(Theta0_b1,w_2(:,5),'b');
scatter(Theta0_b1,Omega_b1,21,'rx')
scatter(Theta0_b2,Omega_b2,21,'gx')
legend('Model 1','Model 2','Experimental Data 1','Experimental Data 2')
xlabel('Angular Position [rad]')
ylabel('Angular Velocity [rad/s]')
title('Angular Velocity vs Angular Position')
hold off


%% Residuals Calculations
% Residual Function Calls
[Res1, Res_mean1, Res_std1, Res1_2, theta1, Res1_mean_2, Res1_std_2] = Residual(Omega_b1,w,Theta0_b1);
[Res2, Res_mean2, Res_std2, Res2_2, theta2, Res2_mean_2, Res2_std_2] = Residual(Omega_b1,w_2(:,5),Theta0_b1);
[Res3, Res_mean3, Res_std3, Res3_2, theta3, Res3_mean_2, Res3_std_2] = Residual(Omega_u2,w_3,Theta0_u2);
[Res4, Res_mean4, Res_std4, Res4_2, theta4, Res4_mean_2, Res4_std_2] = Residual(Omega_u2,w_4,Theta0_u2);

% Standard Deviation and mean calculations
STDs =  [Res1_std_2 Res2_std_2 Res3_std_2 Res4_std_2];
means = [Res1_mean_2 Res2_mean_2 Res3_mean_2 Res4_mean_2];
STD_pos = means + STDs;
STD_neg = means - STDs;

%% Residuals Plots
figure(3)
subplot(2,2,1)
plot(theta1,Res1_2,'k')
hold on
title('Model 1 - Trimmed Residual')
xlabel('Angular Position [rad]')
ylabel('Residual [rad/s]')
ylim([-1.5 .5])
yline(means(1),'--b')
yline(STD_pos(1),'--r')
yline(STD_neg(1),'--r')
legend('Residual', 'Mean Residual','Standard Deviation')
hold off

subplot(2,2,2)
plot(theta2,Res2_2,'k')
hold on
title('Model 2 - Trimmed Residual')
xlabel('Angular Position [rad]')
ylabel('Residual [rad/s]')
ylim([-1 1])
yline(means(2),'--b')
yline(STD_pos(2),'--r')
yline(STD_neg(2),'--r')
legend('Residual', 'Mean Residual','Standard Deviation')
hold off

subplot(2,2,3)
plot(theta3,Res3_2,'k')
hold on
title('Model 3 - Trimmed Residual')
xlabel('Angular Position [rad]')
ylabel('Residual [rad/s]')
ylim([-1 1])
yline(means(3),'--b')
yline(STD_pos(3),'--r')
yline(STD_neg(3),'--r')
legend('Residual', 'Mean Residual','Standard Deviation')
hold off


subplot(2,2,4)
plot(theta4,Res4_2,'k')
hold on
title('Model 4 - Trimmed Residual')
xlabel('Angular Position [rad]')
ylabel('Residual [rad/s]')
ylim([-1 1])
yline(means(4),'--b')
yline(STD_pos(4),'--r')
yline(STD_neg(4),'--r')
legend('Residual', 'Mean Residual','Standard Deviation')
hold off

figure(4)
hold on
plot(Theta0_b1,w,'g');
plot(Theta0_b1,w_2(:,1),'kd-','MarkerSize',3);
plot(Theta0_b1,w_2(:,2),'ko-','MarkerSize',3);
plot(Theta0_b1,w_2(:,3),'ks-','MarkerSize',3);
plot(Theta0_b1,w_2(:,4),'kx-','MarkerSize',3);
plot(Theta0_b1,w_2(:,5),'b*-','MarkerSize',3);
scatter(Theta0_b1,Omega_b1,22,'rx')
legend('Model 1','T = 1 [N*m]','T = 1.2 [N*m]','T = 1.4 [N*m]','T = 1.45 [N*m]','Model 2 : T = 1.5 [N*m] ','Experimental Data')
xlabel('Angular Position [rad]')
ylabel('Angular Velocity [rad/s]')
title('Angular Velocity vs Angular Position')
hold off







