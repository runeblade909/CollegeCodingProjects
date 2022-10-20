function [a,e,i,Omega,omega,f,p] = OE(r,r_dot)

M = 5.97219e24;
G = 6.67408e-11;
mu = G*M;
R = 6371*10^3;

h = cross(r,r_dot);

n = cross([0 0 1],h);


ee =   cross(r_dot,h);
evec = (1/mu).*((ee)-mu.*(r./norm(r)));
e = norm(evec);


% evec = ((norm(r_dot)^2-mu/norm(r))*r - dot(r,r_dot)*r_dot)/mu;
% e = norm(evec);

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

end

