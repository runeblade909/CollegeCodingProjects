%% Zak Reichenbach
%% HW 7

clc
clear all

syms x y z

%Verify Im not insane.....
%3-2-1
M3 = [cos(z) sin(z) 0;...
      -sin(z) cos(z) 0;...
      0 0 1]; 
M2 = [cos(y) 0 -sin(y);...
      0 1 0;...
      sin(y) 0 cos(y)];
M1 = [1 0 0;...
      0 cos(x) sin(x);...
      0 -sin(x) cos(x)];
  
M = M1*M2*M3;


B0 = 0.328474; 
B1 = -0.437966;
B2 = 0.801059;
B3 = -0.242062;

SN = [B0^2+B1^2-B2^2-B3^2, 2*(B1*B2+B0*B3), 2*(B1*B3-B0*B2);...
    2*(B1*B2-B0*B3), B0^2-B1^2+B2^2-B3^2, 2*(B2*B3+B0*B1);...
    2*(B1*B3+B0*B2), 2*(B2*B3-B0*B1), B0^2-B1^2-B2^2+B3^2];

PSI = 230; %deg  SUPPOSED TO BE PSI
theta = 70; %deg
PHI = 103; %dedg

BN = [cosd(theta)*cosd(PSI), cosd(theta)*sind(PSI), -sind(theta);...
      sind(PHI)*sind(theta)*cosd(PSI)-cosd(PHI)*sind(PSI), sind(PHI)*sind(theta)*sind(PHI)+cosd(PHI)*cosd(PSI), sind(PHI)*cosd(theta);...
        cosd(PHI)*sind(theta)*cosd(PSI)+sind(PHI)*sind(PSI), cosd(PHI)*sind(theta)*sind(PHI)-sind(PHI)*cosd(PSI), cosd(PHI)*cosd(theta)];
    
BS = BN*SN';

%Put into E frame with Gamma angle
num = BS(1,1) + BS(2,2) + BS(3,3) + 1;
b0 = 1/2*sqrt(num);
Gamma = 2*acos(b0);

b1 = (BS(2,3)-BS(3,2))/4*b0;

b2 = (BS(3,1)-BS(1,3))/4*b0;

b3 = (BS(1,2)-BS(2,1))/4*b0;


%% Problem 2

%3-1-3 rotation
% (-30, 40, 20)

x = -30;
M31 = [cos(x) sin(x) 0;...
      -sin(x) cos(x) 0;...
      0 0 1]; 
y = 40;
M12 = [1 0 0;...
      0 cos(y) sin(y);...
      0 -sin(y) cos(y)];
z = 20;  
M32= [cos(z) sin(z) 0;...
      -sin(z) cos(z) 0;...
      0 0 1]; 
  
R = M32*M12*M31; 

num2 = R(1,1) + R(2,2) + R(3,3) + 1;
b02 = 1/2*sqrt(num2);
Gamma2 = 2*acos(b02);
Gamma2Prime = Gamma2-2*pi;

mat = [R(2,3)-R(3,2);R(3,1)-R(1,3);R(1,2)-R(2,1)];


E = 1/(2*sin(Gamma2))*mat;
%Check that E is a unit vector;
unit = norm(E);

b12 = (R(2,3)-R(3,2))/4*b02;

b22 = (R(3,1)-R(1,3))/4*b02;

b32 = (R(1,2)-R(2,1))/4*b02;


%% Problem 3
%Body Frame
SB = [0.8273 0.5541 -0.0920]';
MagB = [-0.8285 0.5522 -0.0955]';
%Check
SBCheck = norm(SB);
MagBCheck = norm(MagB);

%Intertial Frame
SN = [-0.1517 -0.9669 0.2050]';
MagN = [-0.8393 0.4494 -0.3044]';

Bt1 = SB;
Bt2 = cross(SB,MagB)/norm(cross(SB,MagB));
Bt3 = cross(Bt1,Bt2);

BT = [Bt1 Bt2 Bt3];

Nt1 = SN;
Nt2 = cross(SN,MagN)/norm(cross(SN,MagN));
Nt3 = cross(Nt1,Nt2);

NT = [Nt1 Nt2 Nt3];

BN = BT*NT';

%3-2-1 Euler angles

psi3 = atand(BN(1,2)/BN(1,1));
theta3 = -asind(BN(1,3));
phi3 = atand(BN(2,3)/BN(3,3));

    
    
    