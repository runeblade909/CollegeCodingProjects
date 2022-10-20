function [rP,vP,rT,vT] = OEPT(r,r_dot,gamma,phi)

M = 5.97219e24;
G = 6.67408e-11;
mu = G*M;
R = 6371*10^3;

h = cross(r,r_dot);

n = cross([0 0 1],h);

evec = ((norm(r_dot)^2-mu/norm(r))*r - dot(r,r_dot)*r_dot)/mu;
e = norm(evec);

energy = norm(r_dot)^2/2-mu/norm(r);

if abs(e-1.0)>eps
   a = -mu/(2*energy);
   p = a*(1-e^2);
else
   p = norm(h)^2/mu;
   a = inf;
end

i = acosd(h(3)/norm(h));

Omega = acosd(n(1)/norm(n));

if n(2)<0
   Omega = 360-Omega;
end
argp = acosd(dot(n,evec)/(norm(n)*e));

if evec(3)<0
   omega = 360-argp;
else
    omega = argp;
end
f = acosd(dot(evec,r)/(e*norm(r)));

if dot(r,r_dot)<0
   f = 360 - f;
end


%Now the perifocal frame

M1 = [1 0 0; 0 cosd(i) sind(i); 0 -sind(i) cosd(i)];

M3 = [cosd(Omega) sind(Omega) 0; -sind(Omega) cosd(Omega) 0; 0 0 1];
M32 = [cosd(omega) sind(omega) 0; -sind(omega) cosd(omega) 0; 0 0 1];

NP = M32*M1*M3;

rP = r*NP;
vP = r_dot*NP;

%Now the topocentric frame


M3T = [cosd(gamma) sind(gamma) 0; -sind(gamma) cosd(gamma) 0; 0 0 1];
M2T = [cosd(phi) 0 -sind(phi); 0 1 0; sind(phi) 0 cosd(phi)];

NT = M3T*M2T;

rT = r*NT;
vT = r_dot*NT;


end

