clc
clear
syms phi theta gama
% phi = 50;
% theta = 25;
% gama = 70;

% 3-1-3
PhiMatrix = [cosd(phi) sind(phi) 0; -sind(phi) cosd(phi) 0; 0 0 1];
ThetaMatrix = [cosd(theta) 0 -sind(theta); 0 1 0 ; sind(theta) 0 cosd(theta)];
GammaMatrix = [cosd(gama) sind(gama) 0; -sind(gama) cosd(gama) 0; 0 0 1];

Phi2Matrix = [1 0 0; 0 cosd(phi) sind(phi); 0 -sind(phi) cosd(phi)];


Q2 = Phi2Matrix*ThetaMatrix*GammaMatrix;



Q = PhiMatrix*ThetaMatrix*GammaMatrix;



angle = acosd(Q(1,1));
