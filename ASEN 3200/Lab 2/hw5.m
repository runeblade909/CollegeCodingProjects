clc
clear all

x = 20;
A = [1 0 0 ; 0 cosd(x) sind(x); 0 -sind(x) cosd(x)];
y = -10;
B =[cosd(y) 0 -sind(y); 0 1 0; sind(y) 0 cosd(y)];
z = 30;
C = [cosd(z) sind(z) 0; -sind(z) cosd(z) 0; 0 0 1];

Rc = 7000*10^3;

R = A*B*C;

Body  = R*[0;0;Rc];

G = 6.67430*10^-11;
Me = 5.9736*10^24;


Const = (3*G*Me)/Rc^5;

Ix = 400;
Iy = 300;
Iz = 200;
L = Const.*[Body(2)*Body(3)*(Iz-Iy);Body(1)*Body(3)*(Ix-Iz);Body(1)*Body(2)*(Iy-Ix)];



%Problem 3

G = 6.67430*10^-11;
Me = 5.9736*10^24;
Mu = G*Me;


Ir = 400;
Ip = 500;
Iy = 300;

R =1;
Const2 = sqrt(Mu/R);
wf = sqrt(3*(Ir-Iy)/Ip);

