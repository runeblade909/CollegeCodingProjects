clc
clear
close all


Brit = xlsread('ASEN1022_Feb21_Brittle_Data.csv');
Duct = xlsread('ASEN1022_Jan28_Ductile_Data.csv');

%now to calculate the stress and strain AT EVERY DATA POINT
DuctStress = zeros(635,1);
BritStress = zeros(1190,1);
DuctStrain = zeros(635,1);
BritStrain = zeros(1190,1);
DuctModulus = zeros(635,1);
BritModulus = zeros(1190,1);

%Im assuming that the 1/2'' and 1/4'' measurements given are the length and
%width of the area we need to be looking at to calculate the stress.

%so to calculate stress, we are going to take the load that varies over
%time in each data set and divide it by the calculated area
AreaDuctile = .504*.187;
AreaBrittle = .500*.191;

for i = 1:635
    DuctStress(i,1) = Duct(i,2)/AreaDuctile;    
end
for i = 1:885
    BritStress(i,1) = Brit(i,2)/AreaBrittle;    
end
for i = 1:635
   DuctStrain(i,1) =  Duct(i,4);
end
for i = 1:885
   BritStrain(i,1) = Brit(i,4); 
end

%---------Brittle Graph-----------------------------------------------------------------------
p2 = polyfit(BritStrain(1:885),BritStress(1:885),1);
x2 = BritStrain(1:885);
y2 = polyval(p2,x2);
%Youngs modulus for the brittle specimin
offset2 = BritStrain(1189)*.002;
figure, scatter(BritStrain(1:1189,1),BritStress(1:1189,1)) % theres a little slip at the beginning so I started at 6 instead of 1
hold on

plot(x2+offset2,y2-5000) %gotta figure out what the y shift is
hold off
title('Brittle Data')
ylabel('Stress (Pounds per square inch)')
xlabel('Strain')

%---------Ductile Graph-------------------------------------------------------------------
p1 = polyfit(DuctStrain(1:145),DuctStress(1:145),1);

x1 = DuctStrain(1:145);
y1 = polyval(p1,x1);
%This created the Young's Modulus for our ductile specimin
offset1 = DuctStrain(634)*.002; %This is our .2% offset for the line
figure, scatter(DuctStrain(1:634,1),DuctStress(1:634,1)) %to shift the graph i added the minuses
hold on

plot(x1+offset1,y1-5000) %line shifted to align with graph more properly with 0.2% offset for modulus
hold off
title('Ductile Data')
ylabel('Stress (Pounds per square inch)')
xlabel('Strain')

%ductile modulus
for i = 1:145
   DuctModulus(i,1) = DuctStress(i,1)/DuctStrain(i,1);
end

%885 for brittle







