clc
clear all
%Zak Reichenbach 10/17/2020
%Professor Jackson Calorimetry Project

%% Import the Data
Data = readtable('SampleD');
Time = table2array(Data(:,2));
TimeMinute = Time/60;
Temp1 = table2array(Data(:,3));
Temp2 = table2array(Data(:,6));
TempWater = table2array(Data(:,4));
TempRoom = table2array(Data(:,5));
T1 = mean(TempWater);
%Average out the 2 Calorimeter's

TrueTemp = (Temp1 + Temp2)/2;

%Set constants for equations and indexing for some reason
Cc = 0.895; %J/gC
MassOfSample = 88.897; % g
MassOfCalorimeter = 0.51; %kg

%error from averaging temperature values
%error in temp1

T1Bar = mean(Temp1);
Temp1Values = zeros(714,1);
for i = 1:714
Temp1Values(i) = (Temp1(i)-T1Bar)^2; 

end
Temp1Summation = sum(Temp1Values);
ErrorInTemp1 = sqrt((1/(714-1))*Temp1Summation);


%error in temp1

T2Bar = mean(Temp2);
Temp2Values = zeros(714,1);
for i = 1:714
Temp2Values(i) = (Temp2(i)-T2Bar)^2; 

end
Temp2Summation = sum(Temp2Values);
ErrorInTemp2 = sqrt((1/(714-1))*Temp2Summation);

%Averaging Error
TrueTempError = sqrt(ErrorInTemp1^2 +ErrorInTemp2^2);
%TrueTempError = TrueTempError/sqrt(2); I think this is only applicable if
%the average is 1 number and it tells us the error for that one average
%number especially if its not weighted.

i = 1;
z = 1;
%% Regression Line Windows
%Get all the time in the window for the first 10 minutes
    for j = 1:714
    if TimeMinute(i) < 10   
    i = i+1;
    end
    end
    N1 = i-1;
    %Get all the time in the window 1:27-1:50
    [TempMax,X1] = max(TrueTemp);
    TempLast = TrueTemp(j);
    
%% Bottom Regression
BottomTime = TimeMinute(1:N1);
BottomTemp = TrueTemp(1:N1);
%Calculating Parts variables for B and M
BottomSumTime = sum(BottomTime);                %Sum X
BottomSumTimeSquared = sum(BottomTime.^2);      %Sum X^2
BottomSumTemp = sum(BottomTemp);                %Sum Y
BottomSumTimeXTemp = sum(BottomTime.*BottomTemp);%Sum X*Y
%Delta
BottomDelta = (N1)*BottomSumTimeSquared - BottomSumTime^2;
%Calculating B and M
BottomB = (BottomSumTimeSquared*BottomSumTemp -BottomSumTime* BottomSumTimeXTemp)/BottomDelta;
BottomM = (N1*BottomSumTimeXTemp - BottomSumTime*BottomSumTemp)/BottomDelta;
%Linear Regression
y1 = BottomM*TimeMinute + BottomB; % Bottom fit line

%In the document it says that to sample is added around 11 minutes in, but
%I think this is inaccurate because after I found the data coresponding to
%11 minutes, it did not seem to meet the requirement for the values and
%objectives the rest of the document asked. So i chose the time to be
%10.127 minutes into the experiment, this is right before we see the
%temperature rise exponentially out of the range it was previously in while at a stable temperature,
%and therefore would make sense to be chosen as the spot right before the sample was added.

%Making Estimates for desired times
BottomTempAtSample = BottomM*10.127 + BottomB;
index = find(TimeMinute < 10.227 & TimeMinute >10.027,1);
TrueTempAtSample = TrueTemp(index);

%% Top Regression
TopTime = TimeMinute(X1:j);
TopTemp = TrueTemp(X1:j);

N2 = j-(X1-1);
%Calculating components for B and M
TopSumTime = sum(TopTime);                %Sum X
TopSumTimeSquared = sum(TopTime.^2);      %Sum X^2
TopSumTemp = sum(TopTemp);                %Sum Y
TopSumTimeXTemp = sum(TopTime.*TopTemp);%Sum X*Y
%Delta
TopDelta = (N2)*TopSumTimeSquared - TopSumTime^2;
%Calculating B and M
TopB = (TopSumTimeSquared*TopSumTemp -TopSumTime* TopSumTimeXTemp)/TopDelta;
TopM = (N2*TopSumTimeXTemp - TopSumTime*TopSumTemp)/TopDelta;
%Linear Regression
y2 = TopM*TimeMinute + TopB; % Top Fit line
%Estimate time
TopTempAtSample = TopM*10.127 + TopB;

%Average minute 10.127 temp and its associated actual temperature spot

AverageTEMP = (TopTempAtSample + BottomTempAtSample)/2;
%Finding the index and its associated time
THESPOT = find(TrueTemp < AverageTEMP+.1 & TrueTemp > AverageTEMP-.1);
TIMESPOT = TimeMinute(THESPOT);
%% Preallocating our T2 and T0 values
T2 = TopM*TIMESPOT + TopB;
T0 = mean(y1);
%% Uncertainties for Temperature

%For loops hand the summation clause of the uncertainty equation
r = zeros(437,1);
for k = 1:437 %the range at which our bottom regression makes its estimates
    r(k) = (y1(k) - TrueTemp(k))^2;
end
r = sum(r);

%I wanna look at the associate error looking from where we exterpilate our
%answer to 

w = zeros(289,1);
for h = 1:289 %the range at which our top regression makes its estimates 
   w(h) = (y2(h+425) - TrueTemp(h+425))^2;
end
w = sum(w);


A = zeros(425,1);
for l = 1:425 %the range at which our top regression makes its estimates 
   A(l) = (TempWater(l) - T1)^2;
end
A = sum(A);

BottomUncertainty = sqrt((1/(437-2))*r); % Uncertainty with T0
TopUncertainty = sqrt((1/(289-2))*w);    % Uncertainty with T2
BoilingWaterUncertainty = sqrt((1/(425-2))*A); %Uncertainty with T1

%% Specific heat

%specific heat of sample =
%(510*0.895*(TopTempAtSample-initial Calorimeter temp))/88.897*(Boiling Water Temp (initial Sample temp) -TopTempAtSample)

Cs = ((MassOfCalorimeter*1000)*Cc*(T2-T0))/(MassOfSample*(T1-T2));

%% General rule with partial derivatives of the specific heat for T0,T1,T2

%This is for the partial derivatives of the specific heat equation
syms x y z  % x = T2  y = T0  z = T1  

f = (x-y)/(z-x);

ddT0 = diff(f,y);
ddT1 = diff(f,z);
ddT2 = diff(f,x);

%Plugging in the actual values for the partial derivatives, and multiplying
%by (Mc*Cc)/Ms since it is constant
dT0 = ((MassOfCalorimeter*1000)*Cc)/(MassOfSample*(T2-T1));                 
dT1 = (-1*((MassOfCalorimeter*1000)*Cc*(T2-T0)))/(MassOfSample*(T2-T1))^2;
dT2 = (((MassOfCalorimeter*1000)*Cc*(T2-T0))/(MassOfSample*(T2-T1))^2) - dT0;

%General Formula for 
UncertaintyInCs = sqrt((dT0*BottomUncertainty)^2+(dT1*BoilingWaterUncertainty)^2+(dT2*TopUncertainty)^2);

%% Error Range (Dotted red lines on the graph)
y1TOP = BottomM*TimeMinute + (BottomB+BottomUncertainty);
y1BOT = BottomM*TimeMinute + (BottomB-BottomUncertainty);

y2TOP = TopM*TimeMinute + (TopB+TopUncertainty);
y2BOT = TopM*TimeMinute + (TopB-TopUncertainty);


%% Plot
figure(1)
plot(TimeMinute,TrueTemp);
hold on;
title('Calorimeter Averaged Data with Error Lines');
ylabel('Temperature (C)')
xlabel('Time (m)')
plot(TimeMinute(1:450),y1(1:450));
plot(TimeMinute(400:714),y2(400:714));

% Error ranges

plot(TimeMinute(1:450),y1TOP(1:450),'--r')
plot(TimeMinute(1:450),y1BOT(1:450),'--r')

plot(TimeMinute(400:714),y2TOP(400:714),'--r')
plot(TimeMinute(400:714),y2BOT(400:714),'--r')

hold off

figure(2)
hold on
title('Calorimeter Averaged Data with Temperature Estimates');
plot(TimeMinute,TrueTemp);
plot(TimeMinute(index),BottomTempAtSample,'r*') % Label T1 on graph
text(TimeMinute(index),BottomTempAtSample,'\leftarrow T0');
plot(TimeMinute(index),TopTempAtSample,'r*')
plot(TimeMinute(index),TrueTempAtSample,'b*')
xline(TimeMinute(index));
xline(TIMESPOT);
yline(AverageTEMP);
plot(TimeMinute(index),AverageTEMP,'c*');
plot(TIMESPOT,AverageTEMP,'g*');
plot(TIMESPOT,T2,'g*'); %Label T2 on the graph
text(TIMESPOT,T2,'\leftarrow T2');
plot(TimeMinute(1:450),y1(1:450));
plot(TimeMinute(400:714),y2(400:714));
