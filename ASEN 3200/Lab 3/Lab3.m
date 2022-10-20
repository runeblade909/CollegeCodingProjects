%Lab 3 Closed Loop Feedback
%Zak Reichenbach
% 11/4/2021

%% House Keeping

clc
clear all
close all

%% Data Load

data{1,1} = load('Torque5_5sec_1secrest');
data{1,2} = load('Torque75_5sec_1secrest');
data{1,3} = load('Torque10_5sec_1secrest');

data{1,4} = load('80K120K2');
data{1,5} = load('80K123K2');
data{1,6} = load('80K125K2');
data{1,7} = load('80K1IdealK2');
data{1,8} = load('77K1IdealK2');


%% Data Analysis
%Time[ms] Torque[mNm] Angular Velocity [rpm] Current[Amp]

%Convert RPM to rad/s
convVel = 2*pi/60;
convTime = 1/1000;

%torque
x1 =linspace(1,100); 

figure()
grid on
subplot(1,2,1)
hold on

%Torque Const
Tc = 25.5; %Nm/amp

%Torques
for i = 1:3
t(:,i) = data{1,i}(:,4).*Tc;
plot(data{1,i}(:,1)*convTime,t)
%Torque vs Time
end

Torque = (mean(t(104:600,:)));
figure(1)
hold off
title('Time Vs. Torque')
xlim([-1 15])
ylim([0 20])
legend('Torque: 5[mNm] Time: 5[s]','Torque: 7.5[mNm] Time: 5[s]','Torque: 10[mNm] Time: 5[s]')
xlabel('Time[sec]')
ylabel('Torque[mNm]')


subplot(1,2,2)
hold on
for i = 1:3
plot(data{1,i}(:,1)*convTime,data{1,i}(:,3))

%Angular Velocity vs Time
end
% plot(data{1,i}(:,1)*convTime,500*data{1,i}(:,1)-470)

for i = 1:3
p = polyfit(data{1,i}(104:600,1)*convTime,data{1,i}(104:600,3)*convVel,1);

accel(i) = p(1,1); %rad/sec^2

f = polyval(p,x1);
end

accel = [accel(3),accel(2),accel(1)];

%Moments of Inertia
I = (Torque/1000)./accel(1,2); %get torque in right units, divide by rad/s^2
Iavg = mean(I);

%Standard Deviation Calc.
num =(I-Iavg).^2;
num = sum(num);
std = sqrt(num/3);

%Subplot 2
legend('Torque: 5[mNm] Time: 5[s]','Torque: 7.5[mNm] Time: 5[s]','Torque: 10[mNm] Time: 5[s]')

hold off
title('Time Vs. Velocity')

xlabel('Time[sec]')
ylabel('Angular Velocity[RPM]')

hold off

%Convert rpm to rad/s
Omeaga = 4000*convVel;

%Max Angular Momentum
MaxL = Omeaga*Iavg; %[Nms]

%% Reaction Wheel Control Response

%data 4-8;
% time[ms] Ref Pos [rad] Measured Pos [rad] Current[Amp] KP KD KI
data{4}(2:end,1) = data{4}(2:end,1) - data{4}(2,1);
data{5}(2:end,1) = data{5}(2:end,1) - data{5}(2,1);
data{6}(2:end,1) = data{6}(2:end,1) - data{6}(2,1);
data{7}(2:end,1) = data{7}(2:end,1) - data{7}(2,1);
data{8}(2:end,1) = data{8}(2:end,1) - data{8}(2,1);

RWTc = 33.5; %mNm/A

%Torque
RWTorque120 = data{4}(:,4)*RWTc; %80K1 120K2
RWTorque123 = data{5}(:,4)*RWTc; %80K1 123K2
RWTorque125 = data{6}(:,4)*RWTc; %80K1 125K2
RWTorqueIK2 = data{7}(:,4)*RWTc; %80K1 Ideal K2
RWTorque77IK2 = data{8}(:,4)*RWTc; %77K1 Ideal K2


figure(2)

subplot(5,2,1)
plot(data{4}(:,1),RWTorque120(:))
xlim([0 15000])
title('RW Torque K1 = 80 K2 = 120')
xlabel('Time[ms]')
ylabel('Torque[mNm]')

subplot(5,2,2)
plot(data{4}(:,1),data{4}(:,3))
hold on 
plot(data{4}(:,1),data{4}(:,2))
yline(max(data{4}(:,3)))
ylim([-.1 max(data{4}(:,3))+.1])
xlim([0 15000])
title('Position K1 = 80 K2 = 120')
xlabel('Time[ms]')
ylabel('Angular Position[rad]')

subplot(5,2,3)
plot(data{5}(:,1),RWTorque123(:))
xlim([0 15000])
title('RW Torque K1 = 80 K2 = 123')
xlabel('Time[ms]')
ylabel('Torque[mNm]')

subplot(5,2,4)
plot(data{5}(:,1),data{5}(:,3))
hold on 
plot(data{5}(:,1),data{5}(:,2))
yline(max(data{5}(:,3)))
ylim([-.1 max(data{5}(:,3))+.1])
xlim([0 15000])
title('Position K1 = 80 K2 = 123')
xlabel('Time[ms]')
ylabel('Angular Position[rad]')

subplot(5,2,5)
plot(data{6}(:,1),RWTorque125(:))
xlim([0 15000])
title('RW Torque K1 = 80 K2 = 125')
xlabel('Time[ms]')
ylabel('Torque[mNm]')

subplot(5,2,6)
plot(data{6}(:,1),data{6}(:,3))
hold on 
plot(data{6}(:,1),data{6}(:,2))
yline(max(data{6}(:,3)))
ylim([-.1 max(data{6}(:,3))+.1])
xlim([0 15000])
title('Position K1 = 80 K2 = 125')
xlabel('Time[ms]')
ylabel('Angular Position[rad]')

subplot(5,2,7)
plot(data{7}(:,1),RWTorqueIK2(:))
xlim([0 15000])
title('RW Torque K1 = 80 K2 = 26.4')
xlabel('Time[ms]')
ylabel('Torque[mNm]')

subplot(5,2,8)
plot(data{7}(:,1),data{7}(:,3))
hold on 
plot(data{7}(:,1),data{7}(:,2))
yline(max(data{7}(:,3)))
ylim([-.1 max(data{7}(:,3))+.1])
xlim([0 15000])
title('Position K1 = 80 K2 = 26.4')
xlabel('Time[ms]')
ylabel('Angular Position[rad]')

subplot(5,2,9)
plot(data{8}(:,1),RWTorque77IK2(:))
xlim([0 15000])
title('RW Torque K1 = 77 K2 = 26.4')
xlabel('Time[ms]')
ylabel('Torque[mNm]')
hold off

subplot(5,2,10)
plot(data{8}(:,1),data{8}(:,3))
hold on 
plot(data{8}(:,1),data{8}(:,2))
yline(max(data{8}(:,3)))
ylim([-.1 max(data{8}(:,3))+.1])
xlim([0 15000])
title('Position K1 = 77 K2 = 26.4')
xlabel('Time[ms]')
ylabel('Angular Position[rad]')






