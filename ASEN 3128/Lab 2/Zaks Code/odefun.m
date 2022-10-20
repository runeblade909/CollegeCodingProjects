function x = odefun(t,s)
   global vCoef f1 f2 f3 f4 g m Ix Iy Iz mu Lc Mc Nc
   
   %S Vector x y z phi theta psi u v w p q r
   
   [X,Y,Z] = aeroDragForce(vCoef,s(7),s(8),s(9));
   [L,M,N] = aeroMoments(s(10),s(11),s(12),mu);
   
   phi = s(4);
   Zc = controlZc(f1,f2,f3,f4,phi,X,Y,Z,m,g);
   
   eq_1 = firstEq(s(4),s(5),s(6),s(7),s(8),s(9));
   eq_2 = secondEq(s(4),s(5),s(6),s(10),s(11),s(12));
   eq_3 = thirdEq(s(4),s(5),s(6),s(10),s(11),s(12),s(7),s(8),s(9),X,Y,Z,Zc,g,m);   
   eq_4 = fourthEq(Ix,Iy,Iz,s(10),s(11),s(12),L,M,N,Lc,Mc,Nc); 
   
   x = [eq_1;eq_2;eq_3;eq_4];
end

function Zc = controlZc(f1,f2,f3,f4,psi,X,Y,Z,m,g)
% %     global  f1 f2 f3 f4 
%     Zc = -f1-f2-f3-f4;
    global m g
    Zc = sqrt(Y^2+(m*g)^2);
    %out = [0;0;-f1-f2-f3-f4];
end

function [X,Y,Z] = aeroDragForce(vCoef,u,v,w)
    X = -vCoef*sqrt(u^2+v^2+w^2)*u;
    Y = -vCoef*sqrt(u^2+v^2+w^2)*v;
    Z = -vCoef*sqrt(u^2+v^2+w^2)*w;
    %out = vCoef*(u^2+v^2+w^2)*[u;v;w];
end

function [L,M,N] = aeroMoments(p,q,r,mu)
global  mu
    L = -mu*(sqrt(p^2+q^2+r^2))*p;
    M = -mu*(sqrt(p^2+q^2+r^2))*q;
    N = -mu*(sqrt(p^2+q^2+r^2))*r;
    %out = -mu*(sqrt(p^2+q^2+r^2))*[p;q;r];    
end

function out = firstEq(phi,theta,psi,uE,vE,wE)
%Linear Velocities
    A = [cos(theta)*cos(psi) sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi) cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi);...
        cos(theta)*sin(psi) sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi) cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi);...
        -sin(theta) sin(phi)*cos(theta) cos(phi)*cos(theta);];
    out = A*[uE;vE;wE];
end

function out = secondEq(phi,theta,psi,p,q,r)
%Angular Rates
    A = [1 sin(phi)*tan(theta) cos(phi)*tan(theta);...
        0 cos(phi) -sin(phi);...
        0 sin(phi)*sec(theta) cos(phi)*sec(theta)];
    out = A*[p;q;r];
end

function out = thirdEq (phi,theta,psi,p,q,r,uE,vE,wE,X,Y,Z,Zc,g,m)
global g m 
%Acceleration
    A = [r*vE - q*wE; p*wE - r*uE; q*uE - p*vE]; %Angular Part
    B = g*[sin(theta); -cos(theta)*sin(phi); -cos(theta)*cos(phi)]; %Gravity Part
    C = [X; Y; Z]/m;                      %Aero Force Part
    D = [0; 0; Zc]/m;   %Control Force Part
    out = A + B + C + D;
%     out(3,1) = 0;
end

function out = fourthEq(Ix,Iy,Iz,p,q,r,L,M,N,Lc,Mc,Nc)
global  Ix Iy Iz  Lc Mc Nc
%Angular Accel. w/ moments
    A = [((Iy-Iz)/Ix)*q*r; ((Iz-Ix)/Iy)*p*r; ((Ix-Iy)/Iz)*p*q];
    B = [L/Ix; M/Iy; N/Iz];
    C = [Lc/Ix; Mc/Iy; Nc/Iz];
    out = A + B + C;
end