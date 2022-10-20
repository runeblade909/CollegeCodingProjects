function [e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N)
% The purpose of this function is to utilized PLLT to output the span
% efficiency factor, Coefficient of lift and Ceoffficient of Induced Drag
% using airfoil characteristics and dimentions.
%
% Input:                             
%                                    
% b = Span                           
% a0_t = Chord Lift Slope            
% a0_r = Root Lift Slope             
% c_t = Tip Chord                    
% c_r = Root Chord                   
% aero_t = Tip Zero-Lift AOA         
% aero_r = Root Zero-Lift AOA        
% geo_t = Tip Geometric AOA          
% geo_r = Root Geometric AOA         
% N = Number of Odd Terms            
%                                    
% Output:                            
%                                    
% e = Span Efficiency Factor         
% c_L = Coefficient of Lift           
% c_Di = Coefficient of Induced Drag 

    % Author: Zak Reichenbach
    % Date: 10/25/2021

%Convert Angles to Radians

aero_t = aero_t*pi/180;
aero_r = aero_r*pi/180; 
geo_t = geo_t*pi/180;
geo_r = geo_r*pi/180;

%Y to Theta Coversion

theta = (1:N)*pi/(2*N);
y = b/2*cos(theta);

%Spanwise Values

a0 = a0_r + (a0_t-a0_r)*y/(b/2);
c = c_r + (c_t-c_r)*y/(b/2);
aero = aero_r + (aero_t-aero_r)*y/(b/2);
geo = geo_r + (geo_t-geo_r)*y/(b/2);

%Evaluate Span and Aspect Ratio

S = 0.5*b*(c_r+c_t);
AR = b^2/S;

%Preallocate the LHS and RHS
A = zeros(N,N);
B = zeros(N,1);


%Fundamental Equation of Prandtl Lifting Line Theory Calculation

for i = 1:N
    for j = 1:N
        A(i,j) = ((4*b)/(a0(i)*c(i)))*sin((2*j-1)*theta(i)) + (2*j-1)*sin((2*j-1)*theta(i))/sin(theta(i));
    end
    B(i) = geo(i) - aero(i);
end

%Solve for Coefficients

Aodd = A\B; %x matrix

% Evaluate Outputs

c_L = Aodd(1)*pi*AR;

delta = 0;
for j = 2:N
    delta = delta + (2*j-1)*(Aodd(j)/Aodd(1))^2;
end
e = 1/(1+delta);

c_Di = (c_L)^2/(pi*e*AR); % c_L^2/(pi*AR)*(1+delta)
