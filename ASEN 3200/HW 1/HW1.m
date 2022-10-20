clc
clear all

x = [.353 .612 -.701; -0.927 .127 .353; .127 -.78 .612];
y = [.951 -.235 -.423; .326 .932 .157; .026 .06 .892 ];

yprime = y';

z = x*yprime;


%5
%syms phi theta psi

phi = 30;
theta = -45;
psi = 60;

phi1 = 10;
theta1 = 25;
psi1 = -15;

T1 = [1 0 0; 0 cosd(phi) sind(phi); 0 -sind(phi) cosd(phi)];
T2= [cosd(theta) 0 -sind(theta);0 1 0; sind(theta) 0 cosd(theta)];
T3 = [cosd(psi) sind(psi) 0; -sind(psi) cosd(psi) 0; 0 0 1];


M1 = [1 0 0; 0 cosd(phi1) sind(phi1); 0 -sind(phi1) cosd(phi1)];
M2= [cosd(theta1) 0 -sind(theta1);0 1 0; sind(theta1) 0 cosd(theta1)];
M3 = [cosd(psi1) sind(psi1) 0; -sind(psi1) cosd(psi1) 0; 0 0 1];

R = T1*T2*T3;
R1 = M1*M2*M3;


Q = R*R1';

BF = [atand(Q(1,2)/Q(1,1)) asind(Q(1,3)) atand(Q(2,3)/Q(3,3))];




