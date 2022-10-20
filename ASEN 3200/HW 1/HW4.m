%House Keeping
clc
clear all
close all

%% Problem 1

%Load Data
coast = load('world_coastline_low.txt');

%Ground tracks

%Planet Stuff

M = 5.97219e24;
G = 6.67408e-11;
mu = G*M;
R = 6371*10^3;


%Given Properties
a = 20000;
e = 0.25;
i = 40;
Omega = 300;
omega = 0;
f = 80;

[s,Lat,Long] = HW4GroundTracks(a,e,i,Omega,omega,f);


for i = 1:length(s)
   Radi(i) = sqrt(s(i,1)^2+s(i,2)^2+s(i,3)^2);   
%    lla(i,:) = eci2lla([s(i,1),s(i,2),s(i,3)],[2000, 1, 1, 17, 17, 17.2]);
end

[Perigee,indexP] = min(Radi);

[Apogee,indexA] = max(Radi);


%Plotting

figure()
plot3(s(:,1),s(:,2),s(:,3))
title('Orbit')
hold on
plot3(0,0,0,'bo')
plot3(s(indexP,1),s(indexP,2),s(indexP,3),'ro')
plot3(s(indexA,1),s(indexA,2),s(indexA,3),'go')
grid on

figure()
title('Ground Track')
plot(coast(:,1),coast(:,2))
hold on
plot(Long,Lat,'o')
grid on
xlim([-180 180])

%% Problem 2


%Inclination Variation 5 TIMES
i = [0,20,40,60,80];

%Given Properties
a = 20000;
e = 0.25;
Omega = 300;
omega = 0;
f = 80;

[sI1,LatI1,LongI1] = HW4GroundTracks(a,e,i(1),Omega,omega,f);
[sI2,LatI2,LongI2] = HW4GroundTracks(a,e,i(2),Omega,omega,f);
[sI3,LatI3,LongI3] = HW4GroundTracks(a,e,i(3),Omega,omega,f);
[sI4,LatI4,LongI4] = HW4GroundTracks(a,e,i(4),Omega,omega,f);
[sI5,LatI5,LongI5] = HW4GroundTracks(a,e,i(5),Omega,omega,f);


for i = 1:length(sI1)
   RadiI1(i) = sqrt(sI1(i,1)^2+sI1(i,2)^2+sI1(i,3)^2);   
end

for i = 1:length(sI2)
   RadiI2(i) = sqrt(sI2(i,1)^2+sI2(i,2)^2+sI2(i,3)^2);  
end

for i = 1:length(sI3)
   RadiI3(i) = sqrt(sI3(i,1)^2+sI3(i,2)^2+sI3(i,3)^2); 
end

for i = 1:length(sI4)
   RadiI4(i) = sqrt(sI4(i,1)^2+sI4(i,2)^2+sI4(i,3)^2); 
end

for i = 1:length(sI5)
   RadiI5(i) = sqrt(sI5(i,1)^2+sI5(i,2)^2+sI5(i,3)^2);
end

[PerigeeI1,indexPI1] = min(RadiI1);
[ApogeeI1,indexAI1] = max(RadiI1);

[PerigeeI2,indexPI2] = min(RadiI2);
[ApogeeI2,indexAI2] = max(RadiI2);

[PerigeeI3,indexPI3] = min(RadiI3);
[ApogeeI3,indexAI3] = max(RadiI3);

[PerigeeI4,indexPI4] = min(RadiI4);
[ApogeeI4,indexAI4] = max(RadiI4);

[PerigeeI5,indexPI5] = min(RadiI5);
[ApogeeI5,indexAI5] = max(RadiI5);


%Plotting

figure()
plot3(sI1(:,1),sI1(:,2),sI1(:,3))
hold on

plot3(sI2(:,1),sI2(:,2),sI2(:,3))
plot3(sI3(:,1),sI3(:,2),sI3(:,3))
plot3(sI4(:,1),sI4(:,2),sI4(:,3))
plot3(sI5(:,1),sI5(:,2),sI5(:,3))

plot3(0,0,0,'bo')

plot3(sI1(indexPI1,1),sI1(indexPI1,2),sI1(indexPI1,3),'ro')
plot3(sI1(indexAI1,1),sI1(indexAI1,2),sI1(indexAI1,3),'go')

plot3(sI2(indexPI2,1),sI2(indexPI2,2),sI2(indexPI2,3),'ro')
plot3(sI2(indexAI2,1),sI2(indexAI2,2),sI2(indexAI2,3),'go')

plot3(sI3(indexPI3,1),sI3(indexPI3,2),sI3(indexPI3,3),'ro')
plot3(sI3(indexAI3,1),sI3(indexAI3,2),sI3(indexAI3,3),'go')

plot3(sI4(indexPI4,1),sI4(indexPI4,2),sI4(indexPI4,3),'ro')
plot3(sI4(indexAI4,1),sI4(indexAI4,2),sI4(indexAI4,3),'go')

plot3(sI5(indexPI5,1),sI5(indexPI5,2),sI5(indexPI5,3),'ro')
plot3(sI5(indexAI5,1),sI5(indexAI5,2),sI5(indexAI5,3),'go')

grid on


figure()
plot(coast(:,1),coast(:,2))
hold on
plot(LongI1,LatI1,'o')
plot(LongI2,LatI2,'*')
plot(LongI3,LatI3,'o')
plot(LongI4,LatI4,'o')
plot(LongI5,LatI5,'o')
grid on
xlim([-180 180])


%Orbital Period Variation 3 TIMES

%Given Properties
a = 20000;
i = 40;
e = 0.25;
Omega = 300;
omega = 0;
f = 80;

[sP1,LatP1,LongP1] = HW4GroundTracksPER(a,e,i,Omega,omega,f,1);
[sP2,LatP2,LongP2] = HW4GroundTracksPER(a,e,i,Omega,omega,f,2);
[sP3,LatP3,LongP3] = HW4GroundTracksPER(a,e,i,Omega,omega,f,3);


for i = 1:length(sP1)
   RadiP1(i) = sqrt(sP1(i,1)^2+sP1(i,2)^2+sP1(i,3)^2);   
end

for i = 1:length(sP2)
   RadiP2(i) = sqrt(sP2(i,1)^2+sP2(i,2)^2+sP2(i,3)^2);  
end

for i = 1:length(sP3)
   RadiP3(i) = sqrt(sP3(i,1)^2+sP3(i,2)^2+sP3(i,3)^2); 
end

[PerigeeP1,indexPP1] = min(RadiP1);
[ApogeeP1,indexAP1] = max(RadiP1);

[PerigeeP2,indexPP2] = min(RadiP2);
[ApogeeP2,indexAP2] = max(RadiP2);

[PerigeeP3,indexPP3] = min(RadiP3);
[ApogeeP3,indexAP3] = max(RadiP3);


%Plotting

figure()
plot3(sP1(:,1),sP1(:,2),sP1(:,3))
hold on

plot3(sP2(:,1),sP2(:,2),sP2(:,3))
plot3(sP3(:,1),sP3(:,2),sP3(:,3))

plot3(0,0,0,'bo')

plot3(sP1(indexPP1,1),sP1(indexPP1,2),sP1(indexPP1,3),'ro')
plot3(sP1(indexAP1,1),sP1(indexAP1,2),sP1(indexAP1,3),'go')

plot3(sP2(indexPP2,1),sP2(indexPP2,2),sP2(indexPP2,3),'ro')
plot3(sP2(indexAP2,1),sP2(indexAP2,2),sP2(indexAP2,3),'go')

plot3(sP3(indexPP3,1),sP3(indexPP3,2),sP3(indexPP3,3),'ro')
plot3(sP3(indexAP3,1),sP3(indexAP3,2),sP3(indexAP3,3),'go')


grid on


figure()
plot(coast(:,1),coast(:,2))
hold on
plot(LongP1,LatP1,'o')
plot(LongP2,LatP2,'*')
plot(LongP3,LatP3,'o')
grid on
xlim([-180 180])


%Eccentricity Variation 5 TIMES
e = [0 .25,.5,.75,.95];

%Given Properties
a = 20000;
i = 40;
Omega = 300;
omega = 0;
f = 80;

[sE1,LatE1,LongE1] = HW4GroundTracks(a,e(1),i,Omega,omega,f);
[sE2,LatE2,LongE2] = HW4GroundTracks(a,e(2),i,Omega,omega,f);
[sE3,LatE3,LongE3] = HW4GroundTracks(a,e(3),i,Omega,omega,f);
[sE4,LatE4,LongE4] = HW4GroundTracks(a,e(4),i,Omega,omega,f);
[sE5,LatE5,LongE5] = HW4GroundTracks(a,e(5),i,Omega,omega,f);


for i = 1:length(s)
   RadiE1(i) = sqrt(sE1(i,1)^2+sE1(i,2)^2+sE1(i,3)^2);   
   RadiE2(i) = sqrt(sE2(i,1)^2+sE2(i,2)^2+sE2(i,3)^2);  
   RadiE3(i) = sqrt(sE3(i,1)^2+sE3(i,2)^2+sE3(i,3)^2);  
   RadiE4(i) = sqrt(sE4(i,1)^2+sE4(i,2)^2+sE4(i,3)^2);  
   RadiE5(i) = sqrt(sE5(i,1)^2+sE5(i,2)^2+sE5(i,3)^2);  
end

[PerigeeE1,indexPE1] = min(RadiE1);
[ApogeeE1,indexAE1] = max(RadiE1);

[PerigeeE2,indexPE2] = min(RadiE2);
[ApogeeE2,indexAE2] = max(RadiE2);

[PerigeeE3,indexPE3] = min(RadiE3);
[ApogeeE3,indexAE3] = max(RadiE3);

[PerigeeE4,indexPE4] = min(RadiE4);
[ApogeeE4,indexAE4] = max(RadiE4);

[PerigeeE5,indexPE5] = min(RadiE5);
[ApogeeE5,indexAE5] = max(RadiE5);


%Plotting

figure()
plot3(sE1(:,1),sE1(:,2),sE1(:,3))
hold on

plot3(sE2(:,1),sE2(:,2),sE2(:,3))
plot3(sE3(:,1),sE3(:,2),sE3(:,3))
plot3(sE4(:,1),sE4(:,2),sE4(:,3))
plot3(sE5(:,1),sE5(:,2),sE5(:,3))

plot3(0,0,0,'bo')

plot3(sE1(indexPE1,1),sE1(indexPE1,2),sE1(indexPE1,3),'ro')
plot3(sE1(indexAE1,1),sE1(indexAE1,2),sE1(indexAE1,3),'go')

plot3(sE2(indexPE2,1),sE2(indexPE2,2),sE2(indexPE2,3),'ro')
plot3(sE2(indexAE2,1),sE2(indexAE2,2),sE2(indexAE2,3),'go')

plot3(sE3(indexPE3,1),sE3(indexPE3,2),sE3(indexPE3,3),'ro')
plot3(sE3(indexAE3,1),sE3(indexAE3,2),sE3(indexAE3,3),'go')

plot3(sE4(indexPE4,1),sE4(indexPE4,2),sE4(indexPE4,3),'ro')
plot3(sE4(indexAE4,1),sE4(indexAE4,2),sE4(indexAE4,3),'go')

plot3(sE5(indexPE5,1),sE5(indexPE5,2),sE5(indexPE5,3),'ro')
plot3(sE5(indexAE5,1),sE5(indexAE5,2),sE5(indexAE5,3),'go')

grid on


figure()
plot(coast(:,1),coast(:,2))
hold on
plot(LongE1,LatE1,'o')
plot(LongE2,LatE2,'*')
plot(LongE3,LatE3,'o')
plot(LongE4,LatE4,'o')
plot(LongE5,LatE5,'o')
grid on
xlim([-180 180])


%% Problem 3




