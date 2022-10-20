%orbital elements calcs from r_ECI, rdot_ECI:
function orbElements = orbitalElements(x, mu)

%vector x in form [r; rdot]; 
r = x(1:3); 
rdot = x(4:6); 

rmag = norm(r);

%specific angular momentum:
h = cross(r,rdot); 
hmag = norm(h); 

%inclination:
i = acos(h(3)/hmag); 

%right ascension of the ascending node:
khat = [0;0;1]; 
N = cross(khat, h);
Nmag = norm(N);

if N(2) >= 0
    
    Omega = acos(N(1)/Nmag); 
else
    Omega = 2*pi - acos(N(1)/Nmag);
end

%eccentricity:
e = 1/mu * (cross(rdot, h) - mu * r./rmag); 
emag = norm(e); 

p = hmag^2/mu;
a = p/(1-emag^2); 

%argument of perigree:
if e(3)>=0
    omega = acos(dot(N./Nmag, e./emag)); 
else
    omega = 2*pi - acos(dot(N./Nmag, e./emag)); 
end

%true anomaly:
vr = dot(rdot,r)./rmag; 

if vr>=0
    theta = acos(dot(e./emag, r./rmag));
else
    theta = 2*pi - acos(dot(e./emag, r./rmag));
end

%storing information, angles in degrees:
orbElements = [a; emag; i; Omega; omega; theta];
end