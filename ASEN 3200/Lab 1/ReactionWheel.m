%Code Poorly Written By: Zak Reichenbach
%9/8/2021
clc
clear all
close all

%Problem 2 Reaction Wheel
data{1,1} = load('2021_XX_XX_TH_UNIT05_RWHEELT#5_t#10');
data{1,2} = load('2021_XX_XX_TH_UNIT05_RWHEELT#6_t#10');
data{1,3} = load('2021_XX_XX_TH_UNIT05_RWHEELT#7_t#10');
data{1,4} = load('2021_XX_XX_TH_UNIT05_RWHEELT#8_t#10');
data{1,5} = load('2021_XX_XX_TH_UNIT05_RWHEELT#10_t#5');
data{1,6} = load('2021_XX_XX_TH_UNIT05_RWHEELT#15_t#3');
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
for i = 1:6
plot(data{1,i}(:,1)*convTime,data{1,i}(:,2))
%Torque vs Time
end
hold off
title('Reaction Wheel Time Vs. Torque')
xlim([-1 15])
ylim([0 20])
legend('Torque: 5[mNm] Time: 10[s]','Torque: 6[mNm] Time: 10[s]','Torque: 7[mNm] Time: 10[s]','Torque: 8[mNm] Time: 10[s]','Torque: 10[mNm] Time: 5[s]','Torque: 15[mNm] Time: 3[s]')
xlabel('Time[sec]')
ylabel('Torque[mNm]')


subplot(1,2,2)
hold on
for i = 1:6
plot(data{1,i}(:,1)*convTime,data{1,i}(:,3))
%Angular Velocity vs Time
end
plot(data{1,i}(:,1)*convTime,500*data{1,i}(:,1)-470)

for i = 1:6
p = polyfit(data{1,i}(100:500,1)*convTime,data{1,i}(100:500,3)*convVel,1);
accel(i) = p(1,1); %rad/sec^2

f = polyval(p,x1);
end

%Torque Const
Tc = 25.5; %Nm/amp

%Torques
t = [5 6 7 8 10 15]; %[mNm] NEED Nm


%Moments of Inertia
I = (t./1000)./accel; %get torque in right units, divide by rad/s^2
Iavg = mean(I([1:3,5,6])); %excluding 4 because its wack

%Standard Deviation Calc.
num =(I-Iavg).^2;
num = sum(num([1:3,5,6]));

std = sqrt(num/5);

legend('Torque: 5[mNm] Time: 10[s]','Torque: 6[mNm] Time: 10[s]','Torque: 7[mNm] Time: 10[s]','Torque: 8[mNm] Time: 10[s]','Torque: 10[mNm] Time: 5[s]','Torque: 15[mNm] Time: 3[s]')

hold off
title('Reaction Wheel Time Vs. Velocity')

xlabel('Time[sec]')
ylabel('Angular Velocity[RPM]')
ylim([-25 6000])
xlim([0 20])
hold off

%Convert rpm to rad/s
Omeaga = 4000*convVel;

%Max Angular Momentum
MaxL = Omeaga*Iavg; %[Nms]
% figure()
% hold on
% for i = 1:6
% p = polyfit(data{1,i}(1:100,1)*convTime,data{1,i}(1:100,3),1);
% f = polyval(p,x1);
% plot(f,x1)
% end
% legend('Torque: 5[mNm] Time: 10[s]','Torque: 6[mNm] Time: 10[s]','Torque: 7[mNm] Time: 10[s]','Torque: 8[mNm] Time: 10[s]','Torque: 10[mNm] Time: 5[s]','Torque: 15[mNm] Time: 3[s]')
% hold off
% 


