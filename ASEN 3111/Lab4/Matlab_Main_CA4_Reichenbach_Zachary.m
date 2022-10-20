% ASEN 3111 - Computational Assignment 4 - Main
% Apply given MATLAB functions to all examples in Chapters 8-10 THAT USE APPENDIX
% A, B, and C TABLES AND FIGURE 9.9 and confirm similar numbers to those obtained
% from tables/chart 9.9 by Anderson.
%Zak Reichenbach
%11/21/21

%House Keeping
clc
clear all 
close all

%A = Isentropic Flow Props
%B = Normal Shock Props
%C = Prandtl-Meyer Func and Mach Angle

%{
% 01 - Ex 8.08 - pg 582
% 02 - Ex 8.11 - pg 594
% 03 - Ex 8.12 - pg 594
% 04 - Ex 8.13 - pg 595
% 05 - Ex 8.14 - pg 597
% 06 - Ex 8.15 - pg 598
% 07 - Ex 8.18 - pg 599
% 08 - Ex 8.21 - pg 601
% 09 - Ex 8.23 - pg 606
% 10 - Ex 9.02 - pg 630
% 11 - Ex 9.03 - pg 631
% 12 - Ex 9.04 - pg 632
% 13 - Ex 9.05 - pg 632
% 14 - Ex 9.06 - pg 636
% 15 - Ex 9.07 - pg 643
% 16 - Ex 9.08 - pg 646
% 17 - Ex 9.09 - pg 654
% 18 - Ex 9.10 - pg 658
% 19 - Ex 9.11 - pg 659
% 20 - Ex 9.12 - pg 663
% 21 - Ex 9.13 - pg 665
% 22 - Ex 10.1 - pg 710
% 23 - Ex 10.2 - pg 711
% 24 - Ex 10.3 - pg 711
% 25 - Ex 10.6 - pg 723
%}

%Constants
R = 287;
Gamma = 1.4;

%Example 8.08
%Calculate p0,T0,T*,a*,M*;
M = 3.5; %M
P = 0.3; %atm
T = 180; %K

[ p0op,t0ot,rho0orho ] =isentropic(M);

p0 = p0op*P;
t0 = t0ot*T;

%for M = 1;
t0otS = 1.2;
Tstar = t0/t0otS;

aStar = sqrt(Gamma*R*Tstar);
a = sqrt(Gamma*R*T);

V = M*a;

Mstar = V/aStar;


fprintf('Example 8.8\n')
fprintf('From Appendix A (Matlab Function), for M = 3.5, p0/p = %f(76.27 App A), T0/T = %f(3.45 App A).\n\n',p0op,t0ot)

%Example 8.11 
%"Consider a normal shock wave in air where the upstream
%flow properties areu1=680 m/s,T1=288 K, andp1=1 atm. 
%Calculate the velocity, temperature, and pressure 
%down-stream of the shock"

u1 = 680; %m/s
T1 = 288; %K

a1 = sqrt(Gamma*R*T1);

M1 = u1/a1;

[ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01 ] = shock_calc( M1 );

fprintf('Example 8.11\n')
fprintf('From Appendix B (Matlab Function), for M1 = %f (2 Calc), p2/p1 = %f (4.5 App B), T2/T1 = %f (1.687 App B).\n\n',M1,p02op01,t2ot1)

% 03 - Ex 8.12 - pg 594
% Normal Shock
%P upstream = 1atm. loss of pressure accross the shock for M1 = 2 M1 = 4
M = 2;
[ p0op,t0ot,rho0orho ] =isentropic(M);

[~,p2op1,~,~,~,p02op01 ] = shock_calc(M);


fprintf('Example 8.12\n')
fprintf('From Appendix A (Matlab Function), for M = %f, p0,1/p1 = %f (7.824 App A), T0/T = %f.\n',M,p0op,t0ot)
fprintf('From Appendix B (Matlab Function), for M = %f, p0,1/p0,2 = %f (0.7209 App B), T0/T = %f.\n\n',M,p02op01,t2ot1)

M = 4;
[ p0op,t0ot,rho0orho ] =isentropic(M);

[~,p2op1,~,~,~,p02op01 ] = shock_calc(M);


fprintf('From Appendix A (Matlab Function), for M = %f, p0,1/p1 = %f (151.8 App A), T0/T = %f.\n',M,p0op,t0ot)
fprintf('From Appendix B (Matlab Function), for M = %f, p0,1/p0,2 = %f (0.1388 App B), T0/T = %f.\n\n',M,p02op01,t2ot1)
% 04 - Ex 8.13 - pg 595

M = 2;
P = 2.65e4;
T = 223.3;
[ p0op,t0ot,rho0orho ] =isentropic(M);
[~,p2op1,~,t2ot1,~,p02op01 ] = shock_calc(M);
fprintf('Example 8.13\n')
fprintf('From Appendix A (Matlab Function), for M = %f, p0/p = %f (7.824 App A), T0/T = %f (1.8 App A).\n',M,p0op,t0ot)
fprintf('From Appendix B (Matlab Function), for M = %f, p0/p = %f (0.7209 App B), T0/T = %f(1.688 App A).\n\n',M,p02op01,t2ot1)

M = 0.2;
P = 1.49e5;
T = 401.9;
[ p0op,t0ot,rho0orho ] =isentropic(M);
[~,p2op1,~,t2ot1,~,p02op01 ] = shock_calc(M);
fprintf('From Appendix A (Matlab Function), for M = %f, p2/p0,2 = %f (1.028 App A), T2/T0,2 = %f (1.008 App A).\n',M,p0op,t0ot)


% 05 - Ex 8.14 - pg 597


M = 10;
P = 2.65e4;
T = 223.3;
[ p0op,t0ot,rho0orho ] =isentropic(M);
[~,p2op1,~,t2ot1,~,p02op01 ] = shock_calc(M);
fprintf('Example 8.14\n')
fprintf('From Appendix A (Matlab Function), for M = %f, p0/p = %f (0.4244x10^5 App A), T0/T = %f (21 App A).\n',M,p0op,t0ot)
fprintf('From Appendix B (Matlab Function), for M = %f, p0/p = %f (0.3045x10^-2 App B), T0/T = %f(1 App A).\n\n',M,p02op01,t2ot1)

M = 0.2;
P = 1.49e5;
T = 401.9;
[ p0op,t0ot,rho0orho ] =isentropic(M);
fprintf('From Appendix A (Matlab Function), for M = %f, p2/p0,2 = %f (1.028 App A), T2/T0,2 = %f (1.008 App A).\n\n',M,p0op,t0ot)

% 06 - Ex 8.15 - pg 598

p2op1 = 4.5;
[ M1 ] = shock_calcGetM1( p2op1 );

[M2,p2op1,rho2orho1,t2ot1,~,p02op01 ] = shock_calc(M1);
fprintf('Example 8.15\n')
fprintf('From Appendix B (Matlab Function), For P2/P1 = 4.5, M = %f (2 App B), M2 = %f (0.5774 App B), rho0/rho = %f (2.667 App B), T0/T = %f(1.687 App B).\n\n',M1,M2,rho2orho1,t2ot1)


% 07 - Ex 8.18 - pg 599

V = 1215; %m/s
T = 300; %K

a = sqrt(Gamma*R*T);
M1 = V/a;

[M2,p2op1,~,t2ot1,~,p02op01 ] = shock_calc(M1);
fprintf('Example 8.18\n')
fprintf('From Appendix B (Matlab Function), For V = 1215 and T = 300, M = %f (3.5 App B), M2 = %f (0.4512 App B), T2/T1 = %f(3.315 App B).\n\n',M1,M2,t2ot1)


% 08 - Ex 8.21 - pg 601

M1 = 3.53;

[M2,p2op1,~,t2ot1,~,p02op01 ] = shock_calc(M1);
fprintf('Example 8.21\n')
fprintf('From Appendix B (Matlab Function), M = %f , M2 = %f (0.4492-Table,0.45-Interpolating,0.45-Calculated).\n\n',M1,M2)


% 09 - Ex 8.23 - pg 606

M1 = 8;

[M2,p2op1,rho2orho1,t2ot1,deltasoR,p02op01  ] = shock_calc(M1);
fprintf('Example 8.23\n')
fprintf('From Appendix B (Matlab Function), M = %f, P2/P1 = %f(82.87 App B), P0,2/P0,1 = %f(8.8488x10^-2).\n\n',M1,p2op1,p02op01)


% 10 - Ex 9.02 - pg 630

M1 = 2;
P = 1;
T = 288;

theta = deg2rad(20); %deg
%Find M, p, T, p0, T0
fprintf('Example 9.2\n')

beta = BTMeq( theta,M1 );
beta = rad2deg(beta);

Mn1 = M1*sind(beta);
[ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01 ] = shock_calc( Mn1 );

fprintf('Beta = %f(53.4 Figure 9.9), From Appendix B (Matlab Function), M = %f (1.606 Calc), P2/P1 = %f(2.82 App B), T2/T1 = %f(1.388 App B).\n',beta,Mn1,p2op1,t2ot1)

[ p0op,t0ot,rho0orho ] =isentropic( M1 );

fprintf('From Appendix A (Matlab Function), M = %f, P2/P1 = %f(7.824 App B), T2/T1 = %f(1.8 App B).\n\n',M1,p0op,t0ot)



% 11 - Ex 9.03 - pg 631


M1 = 2.4;
P = 1;
T = 288;

theta = deg2rad(6.5); %deg
%Find M, p, T, p0, T0
fprintf('Example 9.3\n')

beta = BTMeq(theta,M1);


Mn1 = M1*sin(beta);
[ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01 ] = shock_calc( Mn1 );

fprintf('From Appendix B (Matlab Function), M = %f (1.2 Calc.), P2/P1 = %f(1.513 App B), T2/T1 = %f(1.128 App B).\n\n',Mn1,p2op1,t2ot1)


% 12 - Ex 9.04 - pg 632
p2op1 = 3;
beta = 35;

M1 = shock_calcGetM1(p2op1);
fprintf('Example 9.4\n')
fprintf('From BTM EQ from Figure 9.9, M = %f.\n\n',M1)


% 13 - Ex 9.05 - pg 632

M1 = 3;
beta = 40;
%1) NS 2) OS

[ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01 ] = shock_calc( M1 );

fprintf('Example 9.5\n')
fprintf('From Appendix B (Matlab Function), For case 1 with one Normal Shock, M = %f , P0,2/P0,1 = %f(0.3283 App B).\n',M1,p02op01)

Mn1 = M1*sind(40);

[ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01 ] = shock_calc(Mn1);

fprintf('From Appendix B (Matlab Function), For case 2 with one Oblique shock where beta = 40, M = %f , Mn2 = %f(0.588 App B), P0,2/P0,1 = %f(0.7535 App B).\n',Mn1,M2n,p02op01)

theta = 22;

beta = rad2deg(BTMeq(deg2rad(theta),M1));
fprintf('From BTM EQ from Figure 9.9, with theta = 22, beta = %f.\n',beta)

M2 = M2n/sind(beta-theta);
[ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01 ] = shock_calc(M2);
fprintf('From Appendix B (Matlab Function), For case 2 with one Oblique shock where beta = 40, M = %f (1.9 App B), P0,2/P0,1 = %f(0.7674 App B).\n\n',M2,p02op01)



% 14 - Ex 9.06 - pg 636

%Shock on cone
theta = deg2rad(15);
M = 5;
%Calculate the drag coeff over the wedge
%Drag per unit span Cd = D'/q1S = D'/q1C
%D' = 2p,2Lsin(theta) - 2p,1Lsin(theta)
%L = C/cos(theta)
%Cd = (2ctan(theta))(p2-p1)/q1
%q1 = 1/2rhov^2 = gamma/2*p*M^2

%Cd = (2tan(theta))((p2-p1)/(gamma/2*p1*M^2))



fprintf('Example 9.6\n')
beta = BTMeq(theta,M);
beta = rad2deg(beta);
Mn1 = M*sind(beta);
[M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(Mn1);
fprintf('From BTM EQ from Figure 9.9, with theta = %f, beta = %f.\n',theta,beta)
fprintf('From Appendix B (Matlab Function), M = %f (2.05 Calc) , P2/P1 = %f(4.736 App B).\n\n',Mn1,p2op1)


% 15 - Ex 9.07 - pg 643
%compression corner delfection angle 10
%Mach down stream = 3.6
M = 3.6;
theta = 10;

fprintf('Example 9.7\n')
beta = BTMeq(deg2rad(theta),M);
beta = rad2deg(beta);
Mn1 = M*sind(beta);
[M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(Mn1);
fprintf('From BTM EQ from Figure 9.9, with M1 = %f, theta = %f, beta = %f(20 BTM Chart).\n',M,theta,beta)
fprintf('From Appendix B (Matlab Function), M = %f (1.464 Calc), Mn2 = %f (0.7157 App B), P2/P1 = %f(2.32 App B), T2/T1 = %f(1.294 App B).\n',Mn1,M2n,p2op1,t2ot1)

M2 = M2n/(sind(beta-theta));

%Reflected Shock

beta = BTMeq(deg2rad(theta),M2);
beta = rad2deg(beta);
%This beta is the wave angle (angle between shock and flow)
%Shock angle phi
phi = beta-theta;
Mn1 = M2*sind(beta);
[M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(Mn1);
fprintf('From BTM EQ from Figure 9.9, with M1 = %f, theta = %f, beta = %f(20 BTM Chart).\n',M2,theta,beta)
fprintf('From Appendix B (Matlab Function), M = %f (1.358 Calc), Mn3 = %f (0.7572 App B), P3/P2 = %f(1.991 App B), T3/T2 = %f(1.229 App B).\n\n',Mn1,M2n,p2op1,t2ot1)


% 16 - Ex 9.08 - pg 646
%Detatched curved Bow shock on blunt body
%Wave angle at A = 90 B = 60
%Calculate entorpy behind shock (Straight oblique shock)
BetaA = 90;
BetaB = 60;
Mn1 = 8;

% beta = BTMeq(deg2rad(theta),M2);
% beta = rad2deg(beta);
% %This beta is the wave angle (angle between shock and flow)
% %Shock angle phi
% phi = beta-theta;
% Mn1 = M2*sind(beta);
[M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(Mn1);
% fprintf('From BTM EQ from Figure 9.9, with theta = %f, beta = %f.\n',theta,beta)
fprintf('Example 9.8\n')
fprintf('From Appendix B (Matlab Function), M = %f, P3/P2 = %f(74.5 App B), T3/T2 = %f(13.39 App B).\n',Mn1,p2op1,t2ot1)

Cp = 1004.5;
%This change in entropy is constant along the parabolic lines of the bow
%shock, so you can use this to find the other properties at points on that
%constant entropy bow shock
s2_s1 = Cp*log(t2ot1) - 287*log(p2op1);


Mn1 = Mn1*sind(BetaB);
[M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01] = shock_calc(Mn1);
fprintf('From Appendix B (Matlab Function), M = %f(6.928 Calc), P3/P2 = %f(55.38 App B), T3/T2 = %f(10.2 App B).\n\n',Mn1,p2op1,t2ot1)

% 17 - Ex 9.09 - pg 654

M1 = 1.5;
P1 = 1;
T1 = 288;
%Expansion fan around sharp corner through a deflection angle of 15;
theta = 15;
v1 = PMeq_findv( M1 );
fprintf('Example 9.9\n')
v1 = rad2deg(v1);
fprintf('From Appendix C (Matlab Function), M = 1.5, v1 = %f (11.91 App C).\n',v1)
v2 = deg2rad(v1+theta);
M2 = PMeq_findM(v2,2);
v2 = v1+theta;
fprintf('From Appendix C (Matlab Function), M2 = %f(2.0 App C), v2 = %f (26.91 Calc).\n',M2,v2)


[ p0op,t0ot,rho0orho ] =isentropic(M1);
fprintf('From Appendix A (Matlab Function), M = %f, P0,1/P1 = %f(3.671 App B), T0,1/T1 = %f(1.45 App B).\n',M1,p0op,t0ot)
[ p0op,t0ot,rho0orho ] =isentropic(M2);
fprintf('From Appendix A (Matlab Function), M = %f, P0,1/P1 = %f(7.824 App B), T0,1/T1 = %f(1.8 App B).\n\n',M2,p0op,t0ot)


% 18 - Ex 9.10 - pg 658

%Isentropic Compression surface
M1 = 10;
P1 = 1;
%Compression Turn Angle
theta = 15;
%Calculate Mach and pressure behind wave
v1 = PMeq_findv( M1 );
fprintf('Example 9.10\n')
v1 = rad2deg(v1);
fprintf('From Appendix C (Matlab Function), M = %f, v1 = %f (102.3 App C).\n',M1,v1)
v2 = deg2rad(v1-theta);
M2 = PMeq_findM(v2,2);
v2 = v1-theta;
fprintf('From Appendix C (Matlab Function), M2 = %f(6.4 App C), v2 = %f(87.3 Calc).\n\n',M2,v2)


% 19 - Ex 9.11 - pg 659
%Flow over compression corner
M1 = 10;
P1 = 1;
theta = 15;
%Oblique Shock compression
beta = BTMeq(deg2rad(theta),M1);
beta = rad2deg(beta);
fprintf('Example 9.11\n')
fprintf('From BTM EQ from Figure 9.9, with M1 = %f, theta = %f, beta = %f(20 BTM Chart).\n',M1,theta,beta)
Mn1 = M1*sind(beta);
[ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01 ] = shock_calc( Mn1 );
fprintf('From Appendix B (Matlab Function), M = %f(3.42 Calc),Mn2 = %f, P2/P1 = %f(13.32 App B), P0,2/P0,1 = %f(0.2322 App B).\n',Mn1,M2n,p2op1,p02op01)

M2 = M2n/(sind(beta-theta));

[ p0op,t0ot,rho0orho ] =isentropic( M1 );
fprintf('From Appendix A (Matlab Function), M = %f, P0,1/P1 = %f(4.244x10^4 App A).\n',M1,p0op)

[ p0op,t0ot,rho0orho ] =isentropic( M2 );
fprintf('From Appendix A (Matlab Function), M = %f(5.22 Calc), P0,2/P2 = %f(666.1 App A).\n\n',M2,p0op)


% 20 - Ex 9.12 - pg 663

%Lift and Drag coeff for a flat plat at AOA = 5 in Mach 3
AOA = 5;
M = 3;
%AOA = Theta
v1 = PMeq_findv(M);
v1 = rad2deg(v1);
fprintf('Example 9.12\n')
fprintf('From Appendix C (Matlab Function), M2 = %f, v2 = %f(87.3 App C).\n',M,v1)
v2 = AOA + v1;
M2 = PMeq_findM(deg2rad(v2),2);

[ p0op,t0ot,rho0orho ] =isentropic( M );
fprintf('From Appendix A (Matlab Function), M = %f, P0,2/P2 = %f(36.73 App A).\n',M,p0op)
[ p0op,t0ot,rho0orho ] =isentropic( M2 );
fprintf('From Appendix A (Matlab Function), M = %f (App C), P0,2/P2 = %f(55 App A).\n',M2,p0op)

beta = BTMeq(deg2rad(AOA),M);
beta = rad2deg(beta);
fprintf('From BTM EQ from Figure 9.9, with M1 = %f, theta = %f, beta = %f(23.1 BTM Chart).\n',M,AOA,beta)
Mn1 = M*sind(beta);

[ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01 ] = shock_calc( Mn1 );
fprintf('From Appendix B (Matlab Function), M = %f(1.17 Calc), P3/P1 = %f(1.458 App B).\n\n',Mn1,p2op1)


% 21 - Ex 9.13 - pg 665
%Flatplate and 10 deg included angle wedge at Mach 7.


%A - Flat Plate
fprintf('Example 9.13\n')
fprintf('FLAT PLATE PART A\n')
M1 = 7;
AOA = 10;
v1 = PMeq_findv(M1);
v1 = rad2deg(v1);

fprintf('From Appendix C (Matlab Function), M1 = %f, v2 = %f(90.97 App C).\n',M1,v1)
v2 = AOA + v1;
M2 = PMeq_findM(deg2rad(v2),2);

[ p0op,t0ot,rho0orho ] =isentropic( M1);
fprintf('From Appendix A (Matlab Function), M = %f, P0,2/P2 = %f(.414x10^4 App A).\n',M1,p0op)
[ p0op,t0ot,rho0orho ] =isentropic( M2 );
fprintf('From Appendix A (Matlab Function), M = %f (9.56 App C), P0,2/P2 = %f(0.33x10^5 App A).\n',M2,p0op)

beta = BTMeq(deg2rad(AOA),M1);
beta = rad2deg(beta);
fprintf('From BTM EQ from Figure 9.9, with M1 = %f, theta = %f, beta = %f(16.5 BTM Chart).\n',M1,AOA,beta)
Mn1 = M1*sind(beta);

[ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01 ] = shock_calc( Mn1 );
fprintf('From Appendix B (Matlab Function), Mn1 = %f(2.79 Calc), P3/P1 = %f(4.45 App B),.\n\n',Mn1,p2op1)



%B - Wedge
fprintf('WEDGE PART B\n\n')


%Wave over top of wedge
v2 = v1+5;
M2 = PMeq_findM( deg2rad(v2),2);
fprintf('From Appendix C (Matlab Function), M2 = %f (8.1 App C), v2 = %f (Calcs).\n',M2,v2)
[ p0op,t0ot,rho0orho ] =isentropic( M2 );
fprintf('From Appendix A (Matlab Function), M2 = %f (8.1 App C), P0,2/P2 = %f(1.0897x10^4 App A).\n',M2,p0op)

beta = BTMeq(deg2rad(AOA+5),M1);
beta = rad2deg(beta);
fprintf('From BTM EQ from Figure 9.9, with M1 = %f, theta = %f, beta = %f(23.5 BTM Chart).\n',M1,AOA+5,beta)
Mn1 = M1*sind(beta);

[ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01 ] = shock_calc( Mn1 );
fprintf('From Appendix B (Matlab Function), M = %f(2.79 Calc), P3/P1 = %f(8.915 App B).\n\n',Mn1,p2op1)

% 22 - Ex 10.1 - pg 710
%Isentropic supersonic flow through a convergent-divergent nozzle 
AeTt = 10.25;
P = 5; %atm
T = 600; %R
%Calculate M, rho, T at the nozzle exit
[Msup,Msub] = AoverAstar(1.4,AeTt);
fprintf('Example 10.1 \n')
fprintf('From Appendix A (Matlab Function), M = %f,(3.95 App A)\n\n',Msup)


% 23 - Ex 10.2 - pg 711
%Isentropic Flow Through Convergent-Divergent Nozzle with Exit to Throat
%Ratio of 2.
P = 1;  %atm
T = 288; % K
%Calculate Mach, Presure, Temp, at throuat and exit for 2 cases
%At Throat (Applies to both case b and a) since a is supersonic everywhere
%and b is supersonic in throat

fprintf('Example 10.2 \n')

[ p0op,t0ot,rho0orho ] =isentropic( 1 );
fprintf('From Appendix A (Matlab Function), M = %f, P0,2/P2 = %f(0.528 App A), T0/T = %f(0.833 App A).\n\n',1,1/p0op,1/t0ot)


%a) Flow is supersonic at the exit

[Msup, Msub] = AoverAstar(Gamma,2);
fprintf('From Appendix A (Matlab Function), M = %f,(2.2 App A)\n',Msup)
[ p0op,t0ot,rho0orho ] =isentropic( Msup );
fprintf('From Appendix A (Matlab Function), M = %f, P0,2/P2 = %f(10.69 App A), T0/T = %f(1.968 App A).\n\n',Msup,p0op,t0ot)

%b) Flow is subsonic everywhere but throat, where M=1;
fprintf('From Appendix A (Matlab Function), M = %f,(0.3 App A)\n',Msub)

[ p0op,t0ot,rho0orho ] =isentropic( Msub );
fprintf('From Appendix A (Matlab Function), M = %f, P0,2/P2 = %f(1.064 App A), T0/T = %f(1.018 App A).\n\n',Msub,p0op,t0ot)



% 24 - Ex 10.3 - pg 711
%Nozzle in 10.2 exit pressure is 0.973 atm. Find Mach at Throat and Exit

Pe = 0.973;
Prat = 1/Pe;

[Msub] = isentropicFindM(Gamma,Prat);

fprintf('Example 10.3\n')
fprintf('From Appendix A (Matlab Function), for P0/pe = %f (1.028 Calc), M = %f,(0.2 App A)\n',Prat, Msub)
Arat = AratFindA(Gamma,Msub);

fprintf('From Appendix A (Matlab Function), M = %f (0.2 Calc), Ae/At = %f (2.964 App A).\n',Msub,Arat)
At_As = 1/2 * Arat;

[Msup, Msub] = AoverAstar(Gamma,At_As);
fprintf('From Appendix A (Matlab Function), for At/As = %f (1.482 Calc), M = %f,(0.44 App A)\n',At_As,Msub)



% 25 - Ex 10.6 - pg 723
%Prelim Design of Mach 2 supersonic wind tunnel, calculate ratio of
%diffuser to throat area.

%EQ 10.39 At2/At1 = P01/P02

%Assuming normal shock at diffuser entrance
M1 = 2;
fprintf('Example 10.6\n')
[ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01 ] = shock_calc( M1 );
fprintf('From Appendix B (Matlab Function), M = %f, P0,1/P0,2 = %f(0.7209 App B).\n\n',M1,p02op01)








%% Problem 2
%Do Problem 9.14
%Diamond wedge airfoil from figure 9.36 with half angle of
% e = 10. AOA is 15, calculate lift and wave drag coeff.

%5 sections
%Upstream
%all sides of diamond

%Prandtle Mayer Function Value for Region 2
fprintf('Problem 2: Anderson Problem 9.14 Diamond Wedge Airfoil\n\n')

M = 3;
theta = 5;
v1 = PMeq_findv( M );

fprintf('Calculated Prandtl Mayer Function (App C) value V1 = %f.\n',v1)
v2 = rad2deg(v1) + theta;
fprintf('Calculated Prandtl Mayer Function value V2 = %f.\n',v2)
M2 = PMeq_findM( deg2rad(v2),2 );
fprintf('Mach Number in section 2 (App C) M2 = %f.\n',M2)

[ p0op,t0ot,rho0orho ] =isentropic( M2 );
fprintf('From Appendix A (Matlab Function), M2 = %f (3.27 App C), P0,2/P2 = %f (54.8 App A).\n',M2,p0op)

theta = 20;


v3 = v2 + theta;
fprintf('Calculated Prandtl Mayer Function value V3 = %f.\n',v3)

M3 = PMeq_findM( deg2rad(v3),2 );
fprintf('Mach Number in section 3 (App C) M3 = %f.\n',M3)

[ p0op,t0ot,rho0orho ] =isentropic( M3 );
fprintf('From Appendix A (Matlab Function), M2 = %f (3.27 App C), P0,2/P2 = %f (407.94 App A).\n',M3,p0op)

%Find mach normal to region 4
theta = 25;

beta = BTMeq(deg2rad(theta),M);

beta = rad2deg(beta);
fprintf('From BTM EQ from Figure 9.9, with M1 = %f, theta = %f, beta = %f(44 BTM Chart).\n',M,theta,beta)
Mn1 = M*sind(beta);
fprintf('Shock Normal to Region 4 Mn1 = %f.\n',Mn1)
[ M2n,p2op1,rho2orho1,t2ot1,deltasoR,p02op01 ] = shock_calc( Mn1 );
fprintf('From Appendix B (Matlab Function), M = %f, P0,4/P0,1 = %f(0.6835 App B), P4/P1 = %f (4.881 App B), Mn4 = %f (0.5644 App B).\n',Mn1,p02op01,p2op1,M2n)

M4 = M2n/(sind(beta-theta));
fprintf('Mach Number in section 4 (From Mn2) M4 = %f.\n',M4)

v4 = rad2deg(PMeq_findv(M4));
fprintf('Calculated Prandtl Mayer Function (App C) value V4 = %f.\n',v4)



[ p0op,t0ot,rho0orho ] =isentropic( M4 );
fprintf('From Appendix A (Matlab Function), M4 = %f (1.733 App C), P0,4/P4 = %f (5.189 App A).\n',M4,p0op)

%Prandtl Mayer function value for region 5

theta = 20;

v5 = v4 + theta;

fprintf('Calculated Prandtl Mayer Function value V5 = %f.\n',v5)

M5 = PMeq_findM( deg2rad(v5),2 );
fprintf('Mach Number in section 4 (App C) M4 = %f.\n',M4)

[ p0op,t0ot,rho0orho ] =isentropic( M5 );
fprintf('From Appendix A (Matlab Function), M5 = %f (2.48 App C), P0,5/P5 = %f (16.58 App A).\n',M5,p0op)

[ p0op,t0ot,rho0orho ] =isentropic( M );
fprintf('From Appendix A (Matlab Function), M1 = %f, P0,1/P1 = %f (36.73 App A).\n',M,p0op)

%p2/p1 = p2/p0,2 * p0,2/p0,1 * p0,1/p1

P2oP1 = (1/55.020894)*1*p0op; 

fprintf('P2/P1 = %f\n',P2oP1)

%p3/p1 = p2/p1*p3/p0,3*p0,3/p0,2*p0,2/p2
P3oP1 = P2oP1*(1/406.380519)*1*55.020894;

fprintf('P3/P1 = %f\n',P3oP1)

%Found Earlier from App B
P4oP1 = 4.925008;

fprintf('P4/P1 = %f\n',P4oP1)

%p5/p1 = p5/p0,5*p0,5/p0,4*p0,4/p0,1*p0,1/p1
P5oP1 = 1/16.195589*1*0.679279*36.732722;

fprintf('P5/P1 = %f\n',P5oP1)

% Now to calculate lift and drag FINALLY

fprintf('We know qinf = (Gamma*p1*M1^2)/2\n')

fprintf('and Cl or Cd = (L or D)/(qinf * C)\n')

fprintf('From the figure we know L = p4Lcos25 + p5Lcos5 - p2Lcos5 - p3Lcos25\nWith drag using sin instead of cos.\n')

%From Figure
LoC = 1/(2*cosd(10));

Cl = 2/(Gamma*M^2)*LoC*((P4oP1-P3oP1)*cosd(25) + (P5oP1-P2oP1)*cosd(5));

% D = p4Lsin25 + p5Lsin5 - p2Lsin5 - p3Lsin25

Cd = 2/(Gamma*M^2)*LoC*((P4oP1-P3oP1)*sind(25) + (P5oP1-P2oP1)*sind(5));

fprintf('From combining the previous lines, the Cl = %f, and Cd is %f.\n\n',Cl,Cd)

%% Problem 3
%Recreate Figure 9.9 (BTM Chart)
%Theta on X, Beta on Y, Lines of M;

Mach = linspace(3,10,20);

FA = deg2rad(32);

theta = linspace(0.01,FA,20);

for j = 1:20
    
    for i = 1:20
        beta(j,i) = BTMeq( theta(i),Mach(j));
    end

end

for k = 1:20
plot(rad2deg(theta),rad2deg(beta(k,:))) 
hold on
end

xlabel('Theta [Deg]')
ylabel('Beta [Deg]')
title('Figure 9.9 Beta-Theta-Mach Chart')
