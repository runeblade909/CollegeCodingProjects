function x = propagateState(oe0,t,t_0,MU,J2,Re)
%DESCRIPTION: Computes the propagated position and velocity in km, km/s
%accounting for approximate J2 perturbations
%
%INPUTS:
% oe0       Orbit elements [a,e,i,Om,om,f] at time t0 (km,s,rad)
% t         Current time (s)
% t0        Time at the initial epoch (s)
% MU        Central body's gravitational constant (km^3/s^2)
% J2        Central body's J2 parameter (dimensionless)
% Re        Radius of central body (km)
%
%OUTPUTS:
% x         Position and velocity vectors of the form [r; rdot] (6x1) at
%             time t



%1) Compute the mean orbit elements oe(t) at time t due to J2 perturbations
a = oe0(1);
e = oe0(2);
i = oe0(3);
Om0 = oe0(4);
om0 = oe0(5);

n = sqrt(MU/a^3);

p = a*(1-e^2);

Omega_dot = -3/2*n*J2*(Re/p)^2*cos(i);
omega_dot = 3/2*n*J2*(Re/p)^2*(2-5/2*sin(i)^2);

Omega = Om0 + Omega_dot;
omega = om0 + omega_dot;

M = n*(t-t_0);


%Now the perifocal frame

M1 = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];

M3 = [cos(Omega) sin(Omega) 0; -sin(Omega) cos(Omega) 0; 0 0 1];
M32 = [cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0; 0 0 1];

NP = M32*M1*M3;
PN = NP';


%Newtons Method With Keplers Equation
%2) Solve the time-of-flight problem to compute the true anomaly at time t
error = 1e-8;

if M > pi
    E = M+e/2;
else 
    E = M-e/2;
end

ratio = 1;
while abs(ratio) > error
    ratio = (E - e*sin(E) - M)/(1 - e*cos(E));
    E = E - ratio;
end

     f = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));

%3) Compute r(t), rdot(t) in the perifocal frame
     
r = p/(1+e*cos(f));

rE = r*cos(f);

rP = r*sin(f);

h = sqrt(MU*p);

f_dot = h/r^2;

vE = sqrt(MU/p)*-sin(f);
vP = sqrt(MU/p)*(e+cos(f));
     
%4) Compute r(t), rdot(t) in the ECI frame, save into x
R = PN * [rE;rP;0];
V = PN * [vE;vP;0];
%make sure that function has outputs
x = [R;V];


end








