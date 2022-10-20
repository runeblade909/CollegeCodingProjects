

clear all;
close all;


recuv_tempest;

Va_trim = 22;
h_trim = 2438.5;

wind_inertial = [0;0;0];

trim_definition = [Va_trim; h_trim]; %Position and Velocity


%%% Determine trim
[trim_variables, fval] = CalculateTrimVariables(trim_definition, aircraft_parameters);
%Trim Variables [alpha0; de0; dt0]
[trim_state, trim_input]= TrimStateAndInput(trim_variables, trim_definition);

%% Problem 1a
%%% Linear matrices
[Alon, Blon, Alat, Blat] = AircraftLinearModel(trim_definition, trim_variables, aircraft_parameters);

%% Problem 1b
%Determine Eigenvalues

eigenvaluesLon = eig(Alon);
[eigenvectorsLon,~] = eig(Alon);
%Lambda 1,2 are short period, large magnitude. Lambda 3,4 are phugoid, lower
%magnitude

%find short period and phugoid eigenvectors
spVector = eigenvectorsLon(:,1);
phVector = eigenvectorsLon(:,3);

eigenvaluesLat = eig(Alat);
[eigenvectorsLat,~] = eig(Alat);
%Lambda 2,3 are dutch roll

%calculate natural frequency and damping ratio of each eigenvalue
%{
natFSP_Long = sqrt((-real(eigenvaluesLon(1))).^2 + imag(eigenvaluesLon(1)).^2);
natFPH_Long = sqrt((-real(eigenvaluesLon(3))).^2 + imag(eigenvaluesLon(3)).^2);
dampSP_Long = -real(eigenvaluesLon(1))./natFSP_Long;
dampPH_Long = -real(eigenvaluesLon(3))./natFPH_Long;

natFDR_Lat = sqrt((-real(eigenvaluesLat(2))).^2 + imag(eigenvaluesLat(2)).^2);
dampFDR_Lat = -real(eigenvaluesLat(2))./natFDR_Lat;
%}
fprintf('Longitudinal System Properties');
[wN_Lon,damp_Lon] = damp(Alon);
damp(Alon); %print out

fprintf('Lateral System Properties');
[wN_Lat,damp_Lat] = damp(Alat);
damp(Alat); %print

%not done, find time to half/double
fprintf('Roll time to half\n')
timeCRoll = 0.693/(abs(damp_Lat(1))*wN_Lat(1))
fprintf('Spiral time to double\n')
timeCSpiral = 0.693/(abs(damp_Lat(4))*wN_Lat(4))

%find dutch roll vector

%find full lateral A matrix
AlatFull = [Alat zeros(4,2);0 0 sec(trim_variables(1)) 0 0 0;1 0 0 0 ...
    trim_definition(1)*cos(trim_variables(1)) 0];
[eigvectorsLatFull,~] = eig(AlatFull);

drVector = eigvectorsLatFull(:,4);

%% Problem 1c
%Create phasor plots

%Nondimensionalize x and z velocity
spVector(1) = spVector(1)./trim_definition(1);
spVector(2) = spVector(2)./trim_definition(1);

%Normalize wrt pitch angle
spVector = spVector./spVector(4);

%Graph short period phasor plot
figure()
hold on
plot([0,real(spVector(1))],[0,imag(spVector(1))],'b-o');
plot([0,real(spVector(2))],[0,imag(spVector(2))],'b-o');
plot([0,real(spVector(3))],[0,imag(spVector(3))],'b-o');
plot([0,real(spVector(4))],[0,imag(spVector(4))],'b-o');
title('Short Period Phasor Plot')
xlabel('Real Component')
ylabel('Imaginary Component')
hold off

%Graph phugoid mode
%Nondimensionalize x and z velocity
phVector(1) = phVector(1)./trim_definition(1);
phVector(2) = phVector(2)./trim_definition(1);

%Normalize wrt pitch angle
phVector = phVector./phVector(4);

%Graph short period phasor plot
figure()
hold on
plot([0,real(phVector(1))],[0,imag(phVector(1))],'b-o');
plot([0,real(phVector(2))],[0,imag(phVector(2))],'b-o');
plot([0,real(phVector(3))],[0,imag(phVector(3))],'b-o');
plot([0,real(phVector(4))],[0,imag(phVector(4))],'b-o');
title('Phugoid Phasor Plot')
xlabel('Real Component')
ylabel('Imaginary Component')
hold off

%Lateral Mode Phasor Plot (Dutch Roll)
%scale
drVector = drVector./drVector(4);
figure()
hold on
plot([0,real(drVector(1))],[0,imag(drVector(1))],'b-o');
plot([0,real(drVector(2))],[0,imag(drVector(2))],'b-o');
plot([0,real(drVector(3))],[0,imag(drVector(3))],'b-o');
plot([0,real(drVector(4))],[0,imag(drVector(4))],'b-o');
plot([0,real(drVector(5))],[0,imag(drVector(5))],'b-o');
title('Dutch Roll Phasor Plot')
xlabel('Real Component')
ylabel('Imaginary Component')
hold off

%% Problem 2 Longitudinal Modes
%part a: run phugoid
%%%%% Set initial condition

%phugoid mode first
AlonFull = [Alon zeros(4,2);1 0 0 0 0 0;0 1 0 -trim_definition(1) 0 0];
[AlonV,~] = eig(AlonFull);
%from damping and natural freq, col 3 and 4 are SP, col 5 and 6 are PH

%lon: [u w q theta x_E z_E]'

%get phugoid vector and non-dimensionalize
phVector = real(AlonV(:,5));
scale = deg2rad(4)/phVector(4); %change for 2d
ICLon = phVector*scale;


%Run Linear EOM for Phugoid Mode
sys = ss(AlonFull,[],eye(6),[]);
[yPH,tPH,xPH] = initial(sys,ICLon,250);

%calculate the remaining aircraft variables from initial conditions
X = ones(length(tPH),12);
X(:,1) = X(:,1).*yPH(:,5);
X(:,2) = X(:,2).*trim_state(2);
X(:,3) = X(:,3).*yPH(:,6) - trim_definition(2);
X(:,4) = X(:,4).*trim_state(4);
X(:,5) = X(:,5).*yPH(:,4);
X(:,6) = X(:,6).*trim_state(6);
X(:,7) = X(:,7).*yPH(:,1) + trim_definition(1);
X(:,8) = X(:,8).*trim_state(8);
X(:,9) = X(:,9).*yPH(:,2);
X(:,10) = X(:,10).*trim_state(10);
X(:,11) = X(:,11).*yPH(:,3);
X(:,12) = X(:,12).*trim_state(12);

for i=1:length(tPH)
    uPH(i,:) = trim_input';
end

for i=1:length(tPH)
   wPH(i,:) = wind_inertial'; 
end

%Plot aircraft sim of phugoid mode
%PlotAircraft
PlotAircraftSim(tPH,X,uPH,wPH,'r')

%%% Full sim in ode45
initialC = [trim_state(1) + ICLon(5);trim_state(2);trim_state(3) + ICLon(6);...
    trim_state(4);trim_state(5) + ICLon(4);trim_state(6);...
    trim_state(7) + ICLon(1);trim_state(8);trim_state(9) + ICLon(2);...
    trim_state(10);trim_state(11) + ICLon(3);trim_state(12)];
tfinal = 250;
TSPAN = [0 tfinal];
[TOUT2,YOUT2] = ode45(@(t,y) AircraftEOM(t,y,trim_input,wind_inertial,aircraft_parameters),TSPAN,initialC,[]);


for i=1:length(TOUT2)
    UOUT2(i,:) = trim_input';
end

for i=1:length(TOUT2)
   windOUT(i,:) = wind_inertial'; 
end

%Plot aircraft sim of both modes
PlotAircraftSim(TOUT2,YOUT2,UOUT2,windOUT,'b')

%Part b: run short period
close all;
clear UOUT2 windOUT
%get phugoid vector and non-dimensionalize
spVector = real(AlonV(:,3));
scale = deg2rad(5)/spVector(4); %change for 2d
ICLon = spVector*scale;


%Run Linear EOM for short period Mode
sys = ss(AlonFull,[],eye(6),[]);
[ySP,tSP,xSP] = initial(sys,ICLon,25);

%calculate the remaining aircraft variables from initial conditions
X = ones(length(tSP),12);
X(:,1) = X(:,1).*ySP(:,5);
X(:,2) = X(:,2).*trim_state(2);
X(:,3) = X(:,3).*ySP(:,6) - trim_definition(2);
X(:,4) = X(:,4).*trim_state(4);
X(:,5) = X(:,5).*ySP(:,4);
X(:,6) = X(:,6).*trim_state(6);
X(:,7) = X(:,7).*ySP(:,1) + trim_definition(1);
X(:,8) = X(:,8).*trim_state(8);
X(:,9) = X(:,9).*ySP(:,2);
X(:,10) = X(:,10).*trim_state(10);
X(:,11) = X(:,11).*ySP(:,3);
X(:,12) = X(:,12).*trim_state(12);

for i=1:length(tSP)
    uSP(i,:) = trim_input';
end

for i=1:length(tSP)
   wSP(i,:) = wind_inertial'; 
end

%Plot aircraft sim of short period mode
%PlotAircraft
PlotAircraftSim(tSP,X,uSP,wSP,'r')

%%% Full sim in ode45
initialC = [trim_state(1) + ICLon(5);trim_state(2);trim_state(3) + ICLon(6);...
    trim_state(4);trim_state(5) + ICLon(4);trim_state(6);...
    trim_state(7) + ICLon(1);trim_state(8);trim_state(9) + ICLon(2);...
    trim_state(10);trim_state(11) + ICLon(3);trim_state(12)];
tfinal = 25;
TSPAN = [0 tfinal];
[TOUT2,YOUT2] = ode45(@(t,y) AircraftEOM(t,y,trim_input,wind_inertial,aircraft_parameters),TSPAN,initialC,[]);


for i=1:length(TOUT2)
    UOUT2(i,:) = trim_input';
end

for i=1:length(TOUT2)
   windOUT(i,:) = wind_inertial'; 
end

%Plot aircraft sim of both modes
PlotAircraftSim(TOUT2,YOUT2,UOUT2,windOUT,'b')
% close all;
%% Problem 3: Lateral Mode Plots
%%% Linear simulation
% STUDENTS COMPLETE

%Create Full Lateral Behavior Matrix and find Eigenvectors
AlatFull = [Alat zeros(4,2);0 0 sec(trim_variables(1)) 0 0 0;1 0 0 0 ...
    trim_definition(1)*cos(trim_variables(1)) 0];

[AlatV,~] = eig(AlatFull);

%Compare with ALat for Modes: 4 and 5 are dutch roll, 1 and 2 are, 3 and 6
%are

%Part a: Dutch Roll
%Get dutch roll vector and non-dimensionalize
drVector = real(AlatV(:,5));
scale = deg2rad(2)/drVector(4); %change for 2d
ICLat = drVector*scale;

%run initial response for dutch roll
sys = ss(AlatFull,[],eye(6),[]);
[yDR,tDR,xDR] = initial(sys,ICLat,10);

%calculate aircraft variables from initial conditions and trim state
%guide: [v, p, r, phi, psi, y]
XDR = ones(length(tDR),12);
XDR(:,1) = XDR(:,1).*trim_state(1);
XDR(:,2) = XDR(:,2).*yDR(:,6) + trim_state(2);
XDR(:,3) = XDR(:,3).*trim_state(3);
XDR(:,4) = XDR(:,4).*yDR(:,4) + trim_state(4);
XDR(:,5) = XDR(:,5).*trim_state(5);
XDR(:,6) = XDR(:,6).*yDR(:,5) + trim_state(6);
XDR(:,7) = XDR(:,7).*trim_state(7);
XDR(:,8) = XDR(:,8).*yDR(:,1) + trim_state(8);
XDR(:,9) = XDR(:,9).*trim_state(9);
XDR(:,10) = XDR(:,10).*yDR(:,2) + trim_state(10);
XDR(:,11) = XDR(:,11).*trim_state(11);
XDR(:,12) = XDR(:,12).*yDR(:,3) + trim_state(12);

%create wind vec and control vec
for i=1:length(tDR)
    uDR(i,:) = trim_input';
end

for i=1:length(tDR)
   wDR(i,:) = wind_inertial'; 
end

%Plot dutch roll initial response simulation
PlotAircraftSim(tDR,XDR,uDR,wDR,'r')

%Find response of entire simulation
initialCDR = [trim_state(1);trim_state(2) + ICLat(6);trim_state(3);...
    trim_state(4) + ICLat(4);trim_state(5);trim_state(6) + ICLat(5);...
    trim_state(7);trim_state(8) + ICLat(1);trim_state(9);...
    trim_state(10) + ICLat(2);trim_state(11);trim_state(12) + ICLat(3)];
tfinal = 10;
TSPAN = [0 tfinal];
[TOUT3,YOUT3] = ode45(@(t,y) AircraftEOM(t,y,trim_input,wind_inertial,aircraft_parameters),TSPAN,initialCDR,[]);


for i=1:length(TOUT3)
    UOUT3(i,:) = trim_input';
end

for i=1:length(TOUT3)
   windOUT3(i,:) = wind_inertial'; 
end

%Plot Full simulation of dutch roll
PlotAircraftSim(TOUT3,YOUT3,UOUT3,windOUT3,'b')

% close all;
%Part B Spiral Response
sVector = real(AlatV(:,4));
scale = deg2rad(2)/sVector(4); %change for 2d
ICLat = sVector*scale;

%run initial response for spiral
sys = ss(AlatFull,[],eye(6),[]);
[yS,tS,xS] = initial(sys,ICLat,100);

%calculate aircraft variables from initial conditions and trim state
%guide: [v, p, r, phi, psi, y]
XS = ones(length(tS),12);
XS(:,1) = XS(:,1).*trim_state(1);
XS(:,2) = XS(:,2).*yS(:,6) + trim_state(2);
XS(:,3) = XS(:,3).*trim_state(3);
XS(:,4) = XS(:,4).*yS(:,4) + trim_state(4);
XS(:,5) = XS(:,5).*trim_state(5);
XS(:,6) = XS(:,6).*yS(:,5) + trim_state(6);
XS(:,7) = XS(:,7).*trim_state(7);
XS(:,8) = XS(:,8).*yS(:,1) + trim_state(8);
XS(:,9) = XS(:,9).*trim_state(9);
XS(:,10) = XS(:,10).*yS(:,2) + trim_state(10);
XS(:,11) = XS(:,11).*trim_state(11);
XS(:,12) = XS(:,12).*yS(:,3) + trim_state(12);

%create wind vec and control vec
for i=1:length(tS)
    uS(i,:) = trim_input';
end

for i=1:length(tS)
   wS(i,:) = wind_inertial'; 
end

%Plot spiral initial response simulation
PlotAircraftSim(tS,XS,uS,wS,'r')

%Find response of entire simulation
initialCS = [trim_state(1);trim_state(2) + ICLat(6);trim_state(3);...
    trim_state(4) + ICLat(4);trim_state(5);trim_state(6) + ICLat(5);...
    trim_state(7);trim_state(8) + ICLat(1);trim_state(9);...
    trim_state(10) + ICLat(2);trim_state(11);trim_state(12) + ICLat(3)];
tfinal = 100;
TSPAN = [0 tfinal];
[TOUT4,YOUT4] = ode45(@(t,y) AircraftEOM(t,y,trim_input,wind_inertial,aircraft_parameters),TSPAN,initialCS,[]);


for i=1:length(TOUT4)
    UOUT4(i,:) = trim_input';
end

for i=1:length(TOUT4)
   windOUT4(i,:) = wind_inertial'; 
end

%Plot Full simulation of spiral
PlotAircraftSim(TOUT4,YOUT4,UOUT4,windOUT4,'b')

%Part D: Roll response
% close all;
rVector = real(AlatV(:,6));
scale = deg2rad(2)/rVector(4); %change for 2d
ICLat = rVector*scale;

%run initial response for roll
sys = ss(AlatFull,[],eye(6),[]);
[yR,tR,xR] = initial(sys,ICLat,150);

%calculate aircraft variables from initial conditions and trim state
%guide: [v, p, r, phi, psi, y]
XR = ones(length(tR),12);
XR(:,1) = XR(:,1).*trim_state(1);
XR(:,2) = XR(:,2).*yR(:,6) + trim_state(2);
XR(:,3) = XR(:,3).*trim_state(3);
XR(:,4) = XR(:,4).*yR(:,4) + trim_state(4);
XR(:,5) = XR(:,5).*trim_state(5);
XR(:,6) = XR(:,6).*yR(:,5) + trim_state(6);
XR(:,7) = XR(:,7).*trim_state(7);
XR(:,8) = XR(:,8).*yR(:,1) + trim_state(8);
XR(:,9) = XR(:,9).*trim_state(9);
XR(:,10) = XR(:,10).*yR(:,2) + trim_state(10);
XR(:,11) = XR(:,11).*trim_state(11);
XR(:,12) = XR(:,12).*yR(:,3) + trim_state(12);

%create wind vec and control vec
for i=1:length(tR)
    uR(i,:) = trim_input';
end

for i=1:length(tR)
   wR(i,:) = wind_inertial'; 
end

%Plot roll initial response simulation
PlotAircraftSim(tR,XR,uR,wR,'r')

%Find response of entire simulation
initialCR = [trim_state(1);trim_state(2) + ICLat(6);trim_state(3);...
    trim_state(4) + ICLat(4);trim_state(5);trim_state(6) + ICLat(5);...
    trim_state(7);trim_state(8) + ICLat(1);trim_state(9);...
    trim_state(10) + ICLat(2);trim_state(11);trim_state(12) + ICLat(3)];
tfinal = 150;
TSPAN = [0 tfinal];
[TOUT5,YOUT5] = ode45(@(t,y) AircraftEOM(t,y,trim_input,wind_inertial,aircraft_parameters),TSPAN,initialCR,[]);


for i=1:length(TOUT5)
    UOUT5(i,:) = trim_input';
end

for i=1:length(TOUT5)
   windOUT5(i,:) = wind_inertial'; 
end

%Plot Full simulation of dutch roll
PlotAircraftSim(TOUT5,YOUT5,UOUT5,windOUT5,'b')
