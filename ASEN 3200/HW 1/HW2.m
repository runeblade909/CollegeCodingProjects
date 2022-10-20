%Homework 2
%Zak Reichenbach

%Problem 11.7
clc
clear all



m1 = inertiaMatrix(1,1,1,10);
m2 = inertiaMatrix(-1,-1,-1,10);
m3 = inertiaMatrix(4,-4,4,8);
m4 = inertiaMatrix(-2,2,-2,8);
m5 = inertiaMatrix(3,-3,-3,12);
m6 = inertiaMatrix(-3,3,3,12);

mcg = inertiaMatrix(.267,-.267,.267,60);

Ig = m1+m2+m3+m4+m5+m6-mcg;

%Problem 3

Ia = 1/6*2*2^2;
Ib = 2/5*2;

v = [0 1/sqrt(2) 1/sqrt(2)];
Ic = [18.13 0 0; 0 2.133 0; 0 0 18.13];

Iprime = v*Ic*v';

%Problem 4
I = [70 0 0; 0 70 0; 0 0 15]; %This may be wrong
CP = [1 0 0; 0 .866 .5; 0 -.5 .866]; %This also may be wrong
Ic1 = CP*I*CP';


function mat = inertiaMatrix(x,y,z,mass)

mat = mass* [y^2+z^2 -x*y -x*z; -x*y x^2+z^2 -y*z; -z*x -y*z x^2+y^2];


end



