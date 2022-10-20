%% ASEN 3112 - Lab 3
% Part 2: FEM Results - Resonant Frequencies

% Author: Akanksha Nelacanti & Zak Reichenbach
% Date Created: 04/08/2022
% Last Edited: 04/17/2022

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
freq2 = w2/(2*pi)
freq2 = sort(freq2)
freq2_red = freq2(1:3);

w2_sort = sort(w2);
w2_red = w2_sort(1:3);

%% Four-Element Model
% K_4 * U = (omega^2) * M_4 * U

% Coefficients
c_M4 = rho * A * L / 806400; % [lb-s^2 / in]
c_K4 = 8 * E * I_zz / L^3; % [psi-in]

% Build M4 Matrix 
M4 = c_M4 * [77088, 2916*L,23712,-1284*L, 0, 0, 0, 0, 0, 0;
            2916*L, 172*L^2, 1284*L, -73*L^2, 0, 0, 0, 0, 0, 0;
            23712, 1284*L, 154176, 0, 23712, -1284*L, 0, 0, 0, 0;
            -1284*L, -73*L^2, 0, 344*L^2, 1284*L, -73*L^2, 0, 0, 0, 0;
            0, 0, 23712, 1284*L, 154176, 0, 23712, -1284*L, 0, 0;
            0, 0, -1284*L, -73*L^2, 0, 344*L^2, 1284*L, -73*L^2, 0, 0;
            0, 0, 0, 0, 23712, 1284*L, 154176, 0, 23712, -1284*L;
            0, 0, 0, 0, -1284*L, -73*L^2, 0, 344*L^2, 1284*L, -73*L^2;
            0, 0, 0, 0, 0, 0, 23712, 1284*L, 77088 , -2916*L;
            0, 0, 0, 0, 0, 0, -1284*L, -73*L^2, -2916*L , 172*L^2] ...
          + [0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ;
            0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ;
            0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ;
            0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ;
            0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ;
            0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ;
            0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ;
            0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 ;
            0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , M_T,  S_T;
            0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , S_T,  I_T];
        
K4 = c_K4 * [ 96, 12*L, -96, 12*L, 0, 0, 0, 0, 0, 0;
            12*L, 2*L^2, -12*L,L^2, 0, 0, 0, 0, 0, 0;
            -96, -12*L, 192, 0, -96, 12*L, 0, 0, 0, 0;
            12*L,L^2, 0, 4*L^2, -12*L,L^2, 0, 0, 0, 0;
            0,  0, -96, -12*L, 192, 0, -96, 12*L, 0, 0;
            0, 0, 12*L,L^2, 0, 4*L^2, -12*L,L^2, 0, 0;
            0, 0, 0, 0, -96, -12*L, 192, 0, -96, 12*L;
            0, 0, 0, 0, 12*L,L^2, 0, 4*L^2, -12*L,L^2;
            0, 0, 0, 0, 0, 0, -96, -12*L, 96, -12*L;
            0, 0, 0, 0, 0, 0, 12*L,L^2, -12*L, 2*L^2;];
             
M_red4 = M4(3:end,3:end);
K_red4 = K4(3:end,3:end);

[ev4,evec4]=eig(K_red4,M_red4);
w4 = real(sqrt(diag(evec4)));

% Convert from omega to frequencies [Hz]
freq4 = w4/(2*pi);
freq4 = sort(freq4);
freq4_red = freq4(1:3);

w4_sort = sort(w4);
w4_red = w4_sort(1:3);

%% Print resonant frequencies
fprintf('Eigenvalues for Two Element Method: \n')
fprintf('Eig1 = %0.4f rad/s, Eig2 = %0.4f rad/s, Eig3 = %0.4f rad/s \n\n', w2_red)

fprintf('Resonant Frequencies for Two Element Method: \n')
fprintf('Freq1 = %0.4f Hz, Freq2 = %0.4f Hz, Freq3 = %0.4f Hz \n\n', freq2_red)

fprintf('Eigenvalues for Four Element Method: \n')
fprintf('Eig1 = %0.4f rad/s, Eig2 = %0.4f rad/s, Eig3 = %0.4f rad/s \n\n', w4_red)

fprintf('Resonant Frequencies for Four Element Method: \n')
fprintf('Freq1 = %0.4f Hz, Freq2 = %0.4f Hz, Freq3 = %0.4f Hz \n\n', freq4_red)

%% Question 3 - Zak Reichenbach

%% Provided Code

% % %procedure ploteigenvector (L,ev,ne,nsub,scale);
% % %// declare local variables here if required by language
% % L; 
% % EigenValue2 = [0;0;ev2(:,2)]; 
% % ne=2; 
% % nsub = 10; 
% % scale = 1;
% % nv=ne*nsub+1; 
% % Le=L/ne; 
% % dx=Le/nsub; 
% % k=1;
% % 
% % x = zeros(nv,1);
% % v = zeros(nv,1); %% declare and set to zero plot arrays
% % 
% % for e = 1:ne  %% loop over elements
% %     
% %     xi = Le*(e-1); 
% %     vi = EigenValue2(2*e-1); 
% %     qi = EigenValue2(2*e); 
% %     vj = EigenValue2(2*e+1); 
% %     qj = EigenValue2(2*e+2);
% %     
% %     for n = 1:nsub %% loop over subdivisions
% %         xk = xi+dx*n; 
% %         z = (2*n-nsub)/nsub; %% isoP coordinate
% %         vk = (0.125*(4*(vi+vj)+2*(vi-vj)*(z^2-3)*z + Le*(z^2-1)*(qj-qi+(qi+qj)*z))); %% Hermitian interpolant
% %         k  = k+1;
% %         
% %         x(k)=xk;
% %         v(k)=-vk; %% build plot functions
% %     end  %% end n loop
% %     
% % end %% end e loop
% % 
% % v = v./mean(v);
% % v = scale*v;

%% Plotting response for each frequency eigen mode

[x1,v1] = eigenmode(L,ev2(:,1),2,10,.05);
[x2,v2] = eigenmode(L,ev2(:,2),2,10,1);
[x3,v3] = eigenmode(L,ev2(:,3),2,10,-1);
[x22,v22] = eigenmode(L,ev2(:,4),2,10,1);


ModeFreq1 = freq2(1);
ModeFreq2 = freq2(2);
ModeFreq3 = freq2(3);
ModeFreq22 = freq2(4);

doug = sprintf('Mode Shape 1 %0.1f Hz',ModeFreq1);
dimma = sprintf('Mode Shape 2 %0.1f Hz',ModeFreq2);
dome = sprintf('Mode Shape 3 %0.1f Hz',ModeFreq3);
sone = sprintf('Mode Shape 4 %0.1f Hz',ModeFreq22);

figure
title('Modal Response Two-Element')
hold on
plot(x1,v1)
plot(x2,v2)
plot(x3,v3)
plot(x22,v22)
legend(doug,dimma,dome,sone)
xlim([0 12])
grid on
% yline(0)
% xline(0)
hold off
xlabel('Length(in)')
ylabel('V - Normalized')
%v (vertical) vs x (horizontal) -- language dependent


%% Four-Element Modal Response


[x4,v4] = eigenmode(L,ev4(:,4),4,10,0.05);
[x5,v5] = eigenmode(L,ev4(:,3),4,10,1);
[x6,v6] = eigenmode(L,ev4(:,2),4,10,-1);
[x7,v7] = eigenmode(L,ev4(:,1),4,10,1);
% [x8,v8] = eigenmode(L,ev4(:,5),4,10,1);
% [x9,v9] = eigenmode(L,ev4(:,6),4,10,1);
% [x10,v10] = eigenmode(L,ev4(:,7),4,10,1);
% [x11,v11] = eigenmode(L,ev4(:,8),4,10,1);


ModeFreq4 = freq4(1);
ModeFreq5 = freq4(2);
ModeFreq6 = freq4(3);
ModeFreq7 = freq4(4);
% ModeFreq8 = freq4(5);
% ModeFreq9 = freq4(6);
% ModeFreq10 = freq4(7);
% ModeFreq11 = freq4(8);

figure
hold on

Mode4 = sprintf('Mode Shape 1 %0.1f Hz',ModeFreq4);
Mode5 = sprintf('Mode Shape 2 %0.1f Hz',ModeFreq5);
Mode6 = sprintf('Mode Shape 3 %0.1f Hz',ModeFreq6);
Mode7 = sprintf('Mode Shape 4 %0.1f Hz',ModeFreq7);
% Mode8 = sprintf('Mode Shape 5 %0.1f Hz',ModeFreq8);
% Mode9 = sprintf('Mode Shape 6 %0.1f Hz',ModeFreq9);
% Mode10 = sprintf('Mode Shape 7 %0.1f Hz',ModeFreq10);
% Mode11 = sprintf('Mode Shape 8 %0.1f Hz',ModeFreq11);

plot(x4,v4)
plot(x5,v5)
plot(x6,v6)
plot(x7,v7)
% plot(x8,v8)
% plot(x9,v9)
% plot(x10,v10)
% plot(x11,v11)

legend(Mode4,Mode5,Mode6,Mode7) %,Mode8,Mode9,Mode10,Mode11
xlim([0 12])
title('Modal Response Four-Element')
grid on
hold off
xlabel('Length(in)')
ylabel('V - Normalized')


