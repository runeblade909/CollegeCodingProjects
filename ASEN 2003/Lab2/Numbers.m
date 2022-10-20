%ASEN 2003 Lab 2: Bouncing Balls
%Code by Zak Reichenbach


%% House Keeping
clc
clear
close all

%% Read Data
Trial_1 = table2array(readtable('0756 Data'));
Trial_2 = table2array(readtable('0757 Data'));
Trial_3 = table2array(readtable('0759 Data'));
Trial_4 = table2array(readtable('0760 Data'));
Trial_5 = table2array(readtable('0761 Data'));

%constant of gravity
g=9.8; %[m/s^2]

Range_Height = 1:7;
Range_Time = 1:5;

%% Fix Data
%This is so all the data starts at the same time value, as the Tracker app
%records time from the start of the video, not when data recording starts.

Trial_1(:,1) = Trial_1(:,1)-Trial_1(1,1);
Trial_2(:,1) = Trial_2(:,1)-Trial_2(1,1);
Trial_3(:,1) = Trial_3(:,1)-Trial_3(1,1);
Trial_4(:,1) = Trial_4(:,1)-Trial_4(1,1);
Trial_5(:,1) = Trial_5(:,1)-Trial_5(1,1);

%% Method 1: Height of Bounce
%Find the height of each bounce.

BounceHeight_1 = BounceHeight(Trial_1);
BounceHeight_2 = BounceHeight(Trial_2);
BounceHeight_3 = BounceHeight(Trial_3);
BounceHeight_4 = BounceHeight(Trial_4);
BounceHeight_5 = BounceHeight(Trial_5);

%Hard editing Trials 3 and 4 because god knows why theyre messed up

BounceHeight_1 = [Trial_1(1,3) BounceHeight_1];
BounceHeight_3 = [Trial_3(1,3) 0.52803 BounceHeight_3];
BounceHeight_4 = [Trial_4(1,3) 0.51974 BounceHeight_4];
BounceHeight_5 = [Trial_5(1,3) BounceHeight_5];

%Selecting 7 bounces
BounceHeight_1_T = BounceHeight_1(Range_Height);
BounceHeight_2_T = BounceHeight_2(Range_Height);
BounceHeight_3_T = BounceHeight_3(Range_Height);
BounceHeight_4_T = BounceHeight_4(Range_Height);
BounceHeight_5_T = BounceHeight_5(Range_Height);

%The function E_h provides the Coefficient of restitution from the
%all the bounce heights in a trial
E_Height_1 = E_h(BounceHeight_1_T);
E_Height_2 = E_h(BounceHeight_2_T);
E_Height_3 = E_h(BounceHeight_3_T);
E_Height_4 = E_h(BounceHeight_4_T);
E_Height_5 = E_h(BounceHeight_5_T);

%Average Coefficient of Restitution based on height values
Avg_E_Height = mean([E_Height_1 E_Height_2 E_Height_3 E_Height_4 E_Height_5]);

%% Method 2: Time of Bounces

%Find the time it takes for each bounce
BounceTime_1 = BounceTime(Trial_1);
BounceTime_2 = BounceTime(Trial_2);
BounceTime_3 = BounceTime(Trial_3);
BounceTime_4 = BounceTime(Trial_4);
BounceTime_5 = BounceTime(Trial_5);

%Selecting 5 Bounces
BounceTime_1_T = BounceTime_1(Range_Time);
BounceTime_2_T = BounceTime_2(Range_Time);
BounceTime_3_T = BounceTime_3(Range_Time);
BounceTime_4_T = BounceTime_4(Range_Time);
BounceTime_5_T = BounceTime_5(Range_Time);

%The function E_t provides the Coefficient of restitution from the
%all the bounce times in a trial
E_Time_1 = E_t(BounceTime_1_T);
E_Time_2 = E_t(BounceTime_2_T);
E_Time_3 = E_t(BounceTime_3_T);
E_Time_4 = E_t(BounceTime_4_T);
E_Time_5 = E_t(BounceTime_5_T);

%Average Coefficient of Resitution based on time values
Avg_E_Time = mean([E_Time_1 E_Time_2 E_Time_3 E_Time_4 E_Time_5]);

%% Method 3: Combined Time and Height

Ts_1 = sum(BounceTime_1);
Ts_2 = sum(BounceTime_2);
Ts_3 = sum(BounceTime_3);
Ts_4 = sum(BounceTime_4);
Ts_5 = sum(BounceTime_5);

h0_1 = BounceHeight_1(1);
h0_2 = BounceHeight_2(1);
h0_3 = BounceHeight_3(1);
h0_4 = BounceHeight_4(1);
h0_5 = BounceHeight_5(1);

E_TH_1 = E_th(BounceTime_1,BounceHeight_1);
E_TH_2 = E_th(BounceTime_2,BounceHeight_2);
E_TH_3 = E_th(BounceTime_3,BounceHeight_3);
E_TH_4 = E_th(BounceTime_4,BounceHeight_4);
E_TH_5 = E_th(BounceTime_5,BounceHeight_5);

%Average Coefficient of Resitution based on the combined time and initial height values
Avg_E_TH = mean([E_TH_1 E_TH_2 E_TH_3 E_TH_4 E_TH_5]);


%% Results

fprintf('Our Coefficient of Restitution (e) has be calculated with 3 different models,\n 1st) Height:%f\n 2nd) Time:%f\n 3rd) Time and Height:%f\n',Avg_E_Height,Avg_E_Time,Avg_E_TH)

%% Save Data for Error Propagation

Data_Height = [BounceHeight_1_T' BounceHeight_2_T' BounceHeight_3_T' BounceHeight_4_T' BounceHeight_5_T'];
Data_Time = [BounceTime_1_T' BounceTime_2_T' BounceTime_3_T' BounceTime_4_T' BounceTime_5_T'];

save('Data_Height','Data_Height')
save('Data_Time','Data_Time')

%% Mean and median values for bounce height, time, and initial height
h0 = sort([h0_1 h0_2 h0_3 h0_4 h0_5],'ascend');
Ts = sort([Ts_1 Ts_2 Ts_3 Ts_4 Ts_5],'ascend');
% Heights for each bounce including all trials
h1 = sort([BounceHeight_1_T(2) BounceHeight_2_T(2) BounceHeight_3_T(2) BounceHeight_4_T(2) BounceHeight_5_T(2)],'ascend');
h2 = sort([BounceHeight_1_T(3) BounceHeight_2_T(3) BounceHeight_3_T(3) BounceHeight_4_T(3) BounceHeight_5_T(3)],'ascend');
h3 = sort([BounceHeight_1_T(4) BounceHeight_2_T(4) BounceHeight_3_T(4) BounceHeight_4_T(4) BounceHeight_5_T(4)],'ascend');
h4 = sort([BounceHeight_1_T(5) BounceHeight_2_T(5) BounceHeight_3_T(5) BounceHeight_4_T(5) BounceHeight_5_T(5)],'ascend');
h5 = sort([BounceHeight_1_T(6) BounceHeight_2_T(6) BounceHeight_3_T(6) BounceHeight_4_T(6) BounceHeight_5_T(6)],'ascend');
h6 = sort([BounceHeight_1_T(7) BounceHeight_2_T(7) BounceHeight_3_T(7) BounceHeight_4_T(7) BounceHeight_5_T(7)],'ascend');
% Times for each bounce including all trials
T1 = sort([BounceTime_1_T(1) BounceTime_2_T(1) BounceTime_3_T(1) BounceTime_4_T(1) BounceTime_5_T(1)],'ascend');
T2 = sort([BounceTime_1_T(2) BounceTime_2_T(2) BounceTime_3_T(2) BounceTime_4_T(2) BounceTime_5_T(2)],'ascend');
T3 = sort([BounceTime_1_T(3) BounceTime_2_T(3) BounceTime_3_T(3) BounceTime_4_T(3) BounceTime_5_T(3)],'ascend');
T4 = sort([BounceTime_1_T(4) BounceTime_2_T(4) BounceTime_3_T(4) BounceTime_4_T(4) BounceTime_5_T(4)],'ascend');
T5 = sort([BounceTime_1_T(5) BounceTime_2_T(5) BounceTime_3_T(5) BounceTime_4_T(5) BounceTime_5_T(5)],'ascend');

%Means and Medians
h0_Mean = mean(h0);
h0_Median = median(h0);

Ts_Mean = mean(Ts);
Ts_Median = median(Ts);

h1_Mean = mean(h1);
h1_Median = median(h1);
h2_Mean = mean(h2);
h2_Median = median(h2);
h3_Mean = mean(h3);
h3_Median = median(h3);
h4_Mean = mean(h4);
h4_Median = median(h4);
h5_Mean = mean(h5);
h5_Median = median(h5);

T1_Mean = mean(T1);
T1_Median = median(T1);
T2_Mean = mean(T2);
T2_Median = median(T2);
T3_Mean = mean(T3);
T3_Median = median(T3);
T4_Mean = mean(T4);
T4_Median = median(T4);
T5_Mean = mean(T5);
T5_Median = median(T5);



%% Plot
figure(1)
plot(Trial_1(:,1),Trial_1(:,3))
xlim([0 5])
ylim([0 1])
title('Trial 1')
xlabel('Time [s]')
ylabel('Height [m]')
figure(2)
plot(Trial_2(:,1),Trial_2(:,3))
xlim([0 5])
ylim([0 1])
title('Trial 2')
xlabel('Time [s]')
ylabel('Height [m]')
figure(3)
plot(Trial_3(:,1),Trial_3(:,3))
xlim([0 5])
ylim([0 1])
title('Trial 3')
xlabel('Time [s]')
ylabel('Height [m]')
figure(4)
plot(Trial_4(:,1),Trial_4(:,3))
xlim([0 5])
ylim([0 1])
title('Trial 4')
xlabel('Time [s]')
ylabel('Height [m]')
figure(5)
plot(Trial_5(:,1),Trial_5(:,3))
xlim([0 5])
ylim([0 1])
title('Trial 5')
xlabel('Time [s]')
ylabel('Height [m]')