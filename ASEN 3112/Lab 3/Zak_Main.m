%% ASEN 3112 - Lab 3
% Part 2: FEM Results - Resonant Frequencies

% Author: Akanksha Nelacanti
% Date Created: 04/08/2022
% Last Edited: 04/15/2022

clear
clc
close all
%% Questions

% What units should everything be in?

%% Givens

L = 12; % [in]
L_E = 4.5; % [in]
L_R = 5; % [in]
w = 1; % [in]
h = 1/8; % [in]
h_E = 1/4; % [in]
h_R = 0.040; % [in]
E = 10175000; % [psi]
rho = 0.0002505; % [lb-s^2 / in^4]
M_T = 1.131 * rho; % [lb-s^2 / in]
S_T = 0.5655 * rho; % [lb-s^2]
I_T = 23.124 * rho; % [lb-s^2-in]
A = w * h; % [in^2]
I_zz = (w * h^3) / 12; % [in^4]

%% Two-Element Model
% K_2 * U = (omega^2) * M_2 * U

% Coefficients
c_M2 = (rho * A * L) / 100800; % [lb-s^2 / in]
c_K2 = (4 * E * I_zz) / L^3; % [psi-in]

% Build M2 Matrix
M2 = c_M2 * [19272, 1458*L, 5928, -642*L, 0, 0;...
            1458*L, 172 * (L^2), 642 * L, -73 * (L^2), 0, 0;...
            5928, 642 * L, 38544, 0, 5928, -642 * L;...
            -642 * L, -73 * (L^2), 0, 344 * (L^2), 642 * L, -73 * (L^2);...
            0, 0, -24, -6 * L, 24, -6 * L;...
            0, 0, 6 * L, L^2, -6 * L, 2 * (L^2)];

extra = zeros(6,6);
extra(5,5) = M_T; extra(5,6) = S_T; extra(6,5) = S_T; extra(6,6) = I_T;

M2 = M2 + extra;
M_red2 = M2(3:end, 3:end);

% Build K2 matrix
K2 = c_K2 * [24, 6 * L, -24, 6 * L, 0, 0;...
            6 * L, 2 * (L^2), -6 * L, L^2, 0, 0;...
            -24, -6 * L, 48, 0, -24, 6 * L;...
            6 * L, L^2, 0, 4 * (L^2), -6 * L, L^2;...
            0, 0, -24, -6 * L, 24, -6 * L;...
            0, 0, 6 * L, L^2, -6 * L, 2 * (L^2)];
        
K_red2 = K2(3:end, 3:end);

% Find eigenvalues aka omega
[ev2, evec2] = eig(K_red2, M_red2);
w2 = real(sqrt(diag(evec2)));

% Convert from omega to frequencies [Hz]
freq2 = w2/(2*pi);
freq2 = sort(freq2);
% B = mink(frequencies2, 3 ); % smallest values

%% Four-Element Model
% K_4 * U = (omega^2) * M_4 * U

% Coefficients
c_M4 = (rho * A * L) / 806400; % [lb-s^2 / in]
c_K4 = (8 * E * I_zz) / L^3; % [psi-in]

% Build M4 Matrix

M4 = c_M4 * [
    77088    2916*L    23712    -1284*L       0          0         0          0          0         0
    2916*L   172*L^2   1284*L   -73*L^2       0          0         0          0          0         0
    23712    1284*L    154176      0        23712     -1284*L      0          0          0         0
    -1284*L  -73*L^2      0      344*L^2   1284*L   -73*L^2       0          0          0         0
       0       0       23712     1284*L     154176      0        23712     -1284*L       0         0
       0       0       -1284*L   -73*L^2      0       344*L^2    1284*L    -73*L^2       0         0 
       0       0         0         0        23712     1284*L     154176        0       23712    -1284*L
       0       0         0         0       -1284*L   -73*L^2       0       344*L^2     1284*L   -73*L^2
       0       0         0         0          0          0       23712     1284*L      77088    -2916*L
       0       0         0         0          0          0      -1284*L    -73*L^2    -2916*L    172*L 

];

secoond_part = zeros(10);
secoond_part(9,9)   = M_T;
secoond_part(9,10)  = S_T;
secoond_part(10,9)  = S_T;
secoond_part(10,10) = I_T;

M4 = M4 + secoond_part;

% Build K4 matrix
K4 = c_K4 * [
     96   12*L    -96    12*L     0      0       0        0       0      0
    12*L  2*L^2   -12*L   L^2     0      0       0        0       0      0
    -96  -12*L     192     0     -96    12*L     0        0       0      0
    12*L  L^2       0    4*L^2  -12*L   L*2      0        0       0      0
     0     0       -96   -12*L  192     0       -96     12*L      0      0
     0     0       12*L   L^2    0     4*L^2    -12*L    L^2      0      0
     0     0        0      0    -96    -12*L    192       0      -96    12*L
     0     0        0      0    12*L   L^2       0      4*L^2    -12*L  L^2
     0     0        0      0     0      0       -96     -12*L     96   -12*L
     0     0        0      0     0      0       12*L     L^2    -12*L   2*L
];

M_red4 = M4(3:end, 3:end);
K_red4 = K4(3:end, 3:end);


% Find eigenvalues aka omega
[ev4, evec4] = eig(K_red4, M_red4);
w4 = real(sqrt(diag(evec4)));

% Convert from omega to frequencies [Hz]
freq4 = w4/(2*pi);
freq4 = sort(freq4);

%% My thang

%procedure ploteigenvector (L,ev,ne,nsub,scale);
%// declare local variables here if required by language
L; 
ev2; 
ne=2; 
nsub = 10; 
scale = 1;
nv=ne*nsub+1; 
Le=L/ne; 
dx=Le/nsub; 
k=1;

x = zeros(nv);
v = zeros(nv); %% declare and set to zero plot arrays

for e = 1:ne  %% loop over elements
    
    xi = Le*(e-1); 
    vi = ev2(2*e-1); 
    qi = ev2(2*e); 
    vj = ev2(2*e+1); 
    qj = ev2(2*e+2);
    
    for n = 1:nsub %% loop over subdivisions
        xk = xi+dx*n; 
        z = (2*n-nsub)/nsub; %% isoP coordinate
        vk = scale*(0.125*(4*(vi+vj)+2*(vi-vj)*(z^2-3)*z + Le*(z^2-1)*(qj-qi+(qi+qj)*z))); %% Hermitian interpolant
        k  = k+1;
        
        x(k)=xk;
        v(k)=vk; %% build plot functions
    end  %% end n loop
    
end %% end e loop

figure
plot(x,v)
xlim([0 15])
%v (vertical) vs x (horizontal) -- language dependent





