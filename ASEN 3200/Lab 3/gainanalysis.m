clear
close all
clc

%% load data
datacalc = load('CalculatedGains');
data1 = load('77K1IdealK2');
data2 = load('80K1IdealK2');
data3 = load('80K120K2');
data4 = load('80K123K2');
data5 = load('80K125K2');

%% find start
calcstart = find(diff(datacalc(2:end,2)) > .25,1,'first');
start1 = find(diff(data1(2:end,2)) > .25,1,'first');
start2 = find(diff(data2(2:end,2)) < -.25,1,'first');
start5 = find(diff(data5(2:end,2)) > .25,1,'first');

tcalc = datacalc(calcstart:(calcstart+10/.01),1)/1000;
tcalc = tcalc - min(tcalc);
refcalc = datacalc(calcstart:(calcstart+10/.01),2);
poscalc = datacalc(calcstart:(calcstart+10/.01),3);
currentcalc = datacalc(calcstart:(calcstart+10/.01),4);
Kpcalc = mean(datacalc(calcstart:(calcstart+10/.01),5))/1000;
Kdcalc = mean(datacalc(calcstart:(calcstart+10/.01),6))/1000;
Kicalc = mean(datacalc(calcstart:(calcstart+10/.01),7))/1000;
torquecalc = currentcalc*33.5;

t1 = data1(start1:(start1+10/.01),1)/1000;
t1 = t1 - min(t1);
ref1 = data1(start1:(start1+10/.01),2);
pos1 = data1(start1:(start1+10/.01),3);
current1 = data1(start1:(start1+10/.01),4);
Kp1 = mean(data1(start1:(start1+10/.01),5))/1000;
Kd1 = mean(data1(start1:(start1+10/.01),6))/1000;
Ki1 = mean(data1(start1:(start1+10/.01),7))/1000;
torque1 = current1*33.5;

t2 = data2(start2:(start2+10/.01),1)/1000;
t2 = t2 - min(t2);
ref2 = -1*(data2(start2:(start2+10/.01),2)-.5);
pos2 = -1*(data2(start2:(start2+10/.01),3)-.5);
current2 = data2(start2:(start2+10/.01),4);
Kp2 = mean(data2(start2:(start2+10/.01),5))/1000;
Kd2 = mean(data2(start2:(start2+10/.01),6))/1000;
Ki2 = mean(data2(start2:(start2+10/.01),7))/1000;
torque2 = current2*33.5;

t5 = data5(start5:(start5+10/.01),1)/1000;
t5 = t5 - min(t5);
ref5 = data5(start5:(start5+10/.01),2);
pos5 = data5(start5:(start5+10/.01),3);
current5 = data5(start5:(start5+10/.01),4);
Kp5 = mean(data5(start5:(start5+10/.01),5))/1000;
Kd5 = mean(data5(start5:(start5+10/.01),6))/1000;
Ki5 = mean(data5(start5:(start5+10/.01),7))/1000;
torque5 = current5*33.5;


%% calculate performance
step = mean(diff(tcalc));

% find max over
maxcalc = max(poscalc(1:(10/step)));
maxovercalc = maxcalc - poscalc(round(10/step));
maxperccalc = maxovercalc/.5;

max1 = max(pos1(1:(10/step)));
maxover1 = max1 - pos1(round(10/step));
maxperc1 = maxover1/.5;

max2 = max(pos2(1:(10/step)));
maxover2 = max2 - pos2(round(10/step));
maxperc2 = maxover2/.5;

max5 = max(pos5(1:(10/step)));
maxover5 = max5 - pos5(round(10/step));
maxperc5 = maxover5/.5;

% find settle time
settlecalc1 = tcalc(find((poscalc(1:(10/step)) ...
    > poscalc(round(10/step))*1.05),1,'last')) - tcalc(1);
settlecalc2 = tcalc(find((poscalc(1:(10/step)) ...
    < poscalc(round(10/step))*.95),1,'last')) - tcalc(1);
settlecalc = max([settlecalc1,settlecalc2]);

settle11 = t1(find((pos1(1:(10/step)) ...
    > pos1(round(10/step))*1.05),1,'last')) - t1(1);
settle12 = t1(find((pos1(1:(10/step)) ...
    < pos1(round(10/step))*.95),1,'last')) - t1(1);
settle1 = max([settle11,settle12]);

settle21 = t2(find((pos2(1:(10/step)) ...
    > pos2(round(10/step))*1.05),1,'last')) - t2(1);
settle22 = t2(find((pos2(1:(10/step)) ...
    < pos2(round(10/step))*.95),1,'last')) - t2(1);
settle2 = max([settle21,settle22]);

settle51 = t5(find((pos5(1:(10/step)) ...
    > pos5(round(10/step))*1.05),1,'last')) - t5(1);
settle52 = t5(find((pos5(1:(10/step)) ...
    < pos5(round(10/step))*.95),1,'last')) - t5(1);
settle5 = max([settle51,settle52]);

%% plot

%Simulated Data
figure;
numcalc = Kpcalc/.0066;
dencalc = [1,Kdcalc/.0066,Kpcalc/.0066];
syscalc = tf(numcalc,dencalc);
opt = stepDataOptions('StepAmplitude',.5);
stepplot(syscalc,opt);
hold on
num1 = Kp1/.0066;
den1 = [1,Kd1/.0066,Kp1/.0066];
sys1 = tf(num1,den1);
opt1 = stepDataOptions('StepAmplitude',.5);
stepplot(sys1,opt1);
num2 = Kp2/.0066;
den2 = [1,Kd2/.0066,Kp2/.0066];
sys2 = tf(num2,den2);
opt2 = stepDataOptions('StepAmplitude',.5);
stepplot(sys2,opt2);
num5 = Kp5/.0066;
den5 = [1,Kd5/.0066,Kp5/.0066];
sys5 = tf(num5,den5);
opt5 = stepDataOptions('StepAmplitude',.5);
stepplot(sys5,opt5);
title('Predicted Path');
legend('Calculated Values (K1 = 75.4, K2 = 26.4)'...
    ,'K1 = 77, K2 = Calculated','K1 = 80, K2 = Calculated'...
    ,'K1 = 80, K2 = 25','Location','southeast');
hold off


%Physical Data

figure;

subplot(4,2,1)
hold on
plot(tcalc,poscalc,'b');
plot(tcalc,refcalc,'r');
yline(maxcalc,'k');
xline(settlecalc+tcalc(1),'g');
title('K1 = 75.4 (Calculated Value), K2 = 26.4 (Calculated Value)');
xlabel('Time [s]');
ylabel('Position [rad]');
legend('Actual Path','Reference Height','Maximum Reached','Settling Time'...
    ,'Location','southeast');
axis padded
hold off

subplot(4,2,2)
hold on
plot(tcalc,torquecalc)
title('K1 = 75.4 (Calculated Value), K2 = 26.4 (Calculated Value)');
xlabel('Time [s]');
ylabel('Torque [mNm]');
axis padded
hold off

subplot(4,2,3)
plot(t1,pos1,'b');
hold on
plot(t1,ref1,'r');
yline(max1,'k');
xline(settle1+t1(1),'g');
title('K1 = 77, K2 = 26.4 (Calculated Value)');
xlabel('Time [s]');
ylabel('Position [rad]');
legend('Actual Path','Reference Height','Maximum Reached','Settling Time'...
    ,'Location','southeast');
axis padded
hold off

subplot(4,2,4)
hold on
plot(t1,torque1)
title('K1 = 77, K2 = 26.4 (Calculated Value)');
xlabel('Time [s]');
ylabel('Torque [mNm]');
axis padded
hold off

subplot(4,2,5)
plot(t2,pos2,'b');
hold on
plot(t2,ref2,'r');
yline(max2,'k');
xline(settle2+t2(1),'g');
title('K1 = 80, K2 = 26.4 (Calculated Value)');
xlabel('Time [s]');
ylabel('Position [rad]');
legend('Actual Path','Reference Height','Maximum Reached','Settling Time'...
    ,'Location','southeast');
axis padded
hold off


subplot(4,2,6)
hold on
plot(t2,torque2)
title('K1 = 80, K2 = 26.4 (Calculated Value)');
xlabel('Time [s]');
ylabel('Torque [mNm]');
axis padded
hold off


subplot(4,2,7)
plot(t5,pos5,'b');
hold on
plot(t5,ref5,'r');
yline(max5,'k');
xline(settle5+t5(1),'g');
title('K1 = 80, K2 = 25');
xlabel('Time [s]');
ylabel('Position [rad]');
legend('Actual Path','Reference Height','Maximum Reached','Settling Time'...
    ,'Location','southeast');
axis padded
hold off


subplot(4,2,8)
hold on
plot(t5,torque5)
title('K1 = 80, K2 = 25');
xlabel('Time [s]');
ylabel('Torque [mNm]');
axis padded
hold off

