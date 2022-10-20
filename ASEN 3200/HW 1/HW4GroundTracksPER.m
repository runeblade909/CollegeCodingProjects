function [s,Lat,Long] = HW4GroundTracks(a,e,i,Omega,omega,f,PP)


M = 5.97219e24;
G = 6.67408e-11;
mu = G*M;
R = 6371*10^3;


E = mu/(2*a);   %energy

rp = a*(1-e);   %Radius of Perigee
ra = a*(1+e);   %Radius of Apogee
p = a*(1-e^2);  %Periapse

h = sqrt(p*mu); %Angular Momentum

r = p/(1+e*cos(f)); %Radius

Re = r*cos(f);  %Radius along ie

Rp = r*sin(f);  %Radius along ip

Ve = sqrt(mu/p)*-sin(f);    %Velocity along ie

Vp = sqrt(mu/p)*(e+cos(f)); %Velocity along ip

%we have perifocal frame stuff here
%3-1-3 rotation, Omega, i, omega, then transpose from N->P to P->N
%Earth Centered Inertial

T31 = [cosd(Omega) -sind(Omega) 0; sind(Omega) cosd(Omega) 0; 0 0 1];

T1 = [1 0 0; 0 cosd(i) -sind(i); 0 sind(i) cosd(i)];

T32 = [cosd(omega) -sind(omega) 0; sind(omega) cosd(omega) 0; 0 0 1];

Cnp = T31*T1*T32; %ECI to Parafocal

Cpn = Cnp'; %Parafocal to ECI

rECI = [Re Rp 0]*Cpn;

r_dotECI = [Ve Vp 0]*Cpn;

%Intial Conditions for ODE45
S0 = [rECI r_dotECI];

period = 2*pi*sqrt(a^(3)/mu)*2;

Orbits = 1;

if PP == 1

tspan = [0 period*Orbits];

elseif PP == 2
    
tspan = [0 period*Orbits*2];
   
elseif PP == 3

tspan = [0 period*Orbits*3];
    
end
opts  = odeset('reltol', 1e-12, 'abstol', 1e-12);

[t,s] = ode45(@Satellite,tspan,S0,opts);

wEarth = 7.2921159e-5 * 180/pi; %deg/s

%Convert Cartesian coordinates to spherical ones for Lat and Long

[azimuth,elevation,radius] = cart2sph(s(:,1),s(:,2),s(:,3));

wEarth = (wEarth.*t); %Earth Rotation rate
azimuth = (wEarth.*t) - rad2deg(azimuth);   %Correct Longitude for rotation
elevation = rad2deg(elevation);

%Correct for out of bounds stuff
for i = 1:length(azimuth)
    if azimuth(i) > 180
        azimuth(i) = azimuth(i) - 360;
    end
    if azimuth(i) < -180
       azimuth(i) = azimuth(i) + 360; 
    end
end

Lat = elevation;
Long = azimuth;
end

