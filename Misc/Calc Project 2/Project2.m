%APPM 2350 Project 2
%Zak Reichenbach
%11/9/2019

%Bach Peak         z = e^(-(1.1*sqrt((x-1)^(2)+ (y-3.5)^(2)))^(2))
%W=2000

%Mount Beethoven   z = e^(-(2*sqrt((x-4.5)^(2)+ (y-1.5)^(2)))^(2))
%W=2000

%Hayde Valley      z = e^(-(1*sqrt((x-2.8)^(2)+ (y-3.7)^(2)))^(2))
%W=-950

%Vivaldi Ridge     z = (1)/(1+ (1.5*sqrt((x-2.5)^(2)+ (y-1.8)^(2)))^(2))
%W=1000

%Mount Mozart      z = (1)/(1+ (1.1*sqrt((x-3)^(2)+ (y-2)^(2)))^(2))
%W=1500

%Handel Mountain   z = (1)/(sqrt(1+ (1.9*sqrt((x-1)^(2)+ (y-1)^(2)))^(2)))
%W=3500

%Schumann Hill     z = (1)/(sqrt(1+ (2*sqrt((x-4.25)^(2)+ (y-4.1)^(2)))^(2)))
%W=2000

%Praetorius Valley z = sqrt(1+ (0.5*sqrt((x-2.75)^(2)+ (y-3.75)^(2)))^(2))
%W=200

%Main Plain        z = 1
%W=7000

%Combines Equation z = e^(-(2*sqrt((x-4.5)^(2)+(y-1.5)^(2)))^(2))+e^(-sqrt((x-2.8)^(2)+(y-3.7)^(2))^(2))+(1)/(1+ (1.5*sqrt((x-2.5)^(2)+ (y-1.8)^(2)))^(2))+(1)/(1+ (1.1*sqrt((x-3)^(2)+ (y-2)^(2)))^(2))+(1)/(sqrt(1+ (1.9*sqrt((x-1)^(2)+ (y-1)^(2)))^(2)))+(1)/(sqrt(1+ (2*sqrt((x-4.25)^(2)+ (y-4.1)^(2)))^(2)))+ sqrt(1+ (0.5*sqrt((x-2.75)^(2)+ (y-3.75)^(2)))^(2))+e^(-(1.1*sqrt((x-1)^(2)+ (y-3.5)^(2)))^(2))

%Combined Equations with SCALED weight: z = 2e^(-(2*sqrt((x-4.5)^(2)+ (y-1.5)^(2)))^(2))-0.95e^(-sqrt((x-2.8)^(2)+ (y-3.7)^(2))^(2))+(1)/(1+ (1.5*sqrt((x-2.5)^(2)+ (y-1.8)^(2)))^(2))+(1.5)/(1+ (1.1*sqrt((x-3)^(2)+ (y-2)^(2)))^(2))+(3.5)/(sqrt(1+ (1.9*sqrt((x-1)^(2)+ (y-1)^(2)))^(2)))+(2)/(sqrt(1+ (2*sqrt((x-4.25)^(2)+ (y-4.1)^(2)))^(2)))+0.2*sqrt(1+ (0.5*sqrt((x-2.75)^(2)+ (y-3.75)^(2)))^(2))+2e^(-(1.1*sqrt((x-1)^(2)+ (y-3.5)^(2)))^(2))

clc
clear
[x,y] = meshgrid(0:0.5:5,0:5);
BackPeak = gaussian(1,3.5,1.1);
disp(BackPeak);


% [x,y] = meshgrid(1:0.5:5,1:5);
% z = 2*exp(-(2*sqrt((x-4.5)^(2)+ (y-1.5)^(2)))^(2))-0.95*exp(-sqrt((x-2.8)^(2)+ (y-3.7)^(2))^(2))+(1)/(1+ (1.5*sqrt((x-2.5)^(2)+ (y-1.8)^(2)))^(2))+(1.5)/(1+ (1.1*sqrt((x-3)^(2)+ (y-2)^(2)))^(2))+(3.5)/(sqrt(1+ (1.9*sqrt((x-1)^(2)+ (y-1)^(2)))^(2)))+(2)/(sqrt(1+ (2*sqrt((x-4.25)^(2)+ (y-4.1)^(2)))^(2)))+0.2*sqrt(1+ (0.5*sqrt((x-2.75)^(2)+ (y-3.75)^(2)))^(2))+2*exp(-(1.1*sqrt((x-1)^(2)+ (y-3.5)^(2)))^(2));
% surfc(x,y,z)