%3200 Attitude Dynamics & Control
%Zak Reichenbach
%Pete Dillman
%Alex Lam
%

clc
clear all
close all

          %Initial Speed [km/hr], Procession Time, Final Speed [km/hr]
BikeData = [37 7.35 33;...
            28 4.59 26;...
            50 9.63 44;...
            35 6.38 33];
D = 64; %[cm]
AxelLength = 17; %[cm]
Mass = 2.98; %[kg]
        
%Omega = V/R
%MOI  1)Disk 2)Ring
%Precession rate equation for psi dot with an axial symmetric object
     
%Converting linear speed to meters per second   km/h  (m/km)(h/min)(min/sec)
LinVel = BikeData(:,1)*(1000)*(1/(60*60)); %[m/s]

AngVel = LinVel/(D/2);

MOI_Disc = 1/2*Mass*(D/2)^2;
MOI_Ring = Mass*(D/2)^2;

H_Disc = AngVel*MOI_Disc;
H_Ring = AngVel*MOI_Ring;

psiDot = H_Disc/MOI_Disc;


%L = AxelLength *M*g;

%ThetaDot = L/H