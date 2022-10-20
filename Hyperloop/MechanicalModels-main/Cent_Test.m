%% Housekeeping

clc; clear; close all;

ring_diameter = .5173; % center diameter, m 
circumference = ring_diameter*pi; % center circumference, m
% E_rod = ; % Youngs modulus of steel rod, Pa
A = 2*pi*(ring_diameter/2)^2; % Surface area, m^2

%% Info, constants
rho_sand = 101.82*16.0185; % kg/m^3
rho_tungst = 0.70809488 * 12^3 * 16.0185; % kg/m^3

% Test Tunnel segment dimensions
E_tarp = 0.29*10^9; % Youngs Modulus of tarp, Pa
L_total = 13.5*2.54/100; % Length of whole segment in meters
L_rings = L_total/2; % Separation center-to-center between two rings
R_tunnel = 3*2.54/100; % Radius of a ring, calculated from rod length [m]
Circ_tunnel = 2*pi*R_tunnel; % Circumference of tunnel in meters

%%

%% Modelling force on tunnel
N = 100;
H_INIT = 0.1; % 10 cm from top of tunnel to top of sand

X = linspace(-R_tunnel, R_tunnel, N);
for i = 1:2*N-1
    X = [X linspace(-R_tunnel, R_tunnel, N)];
    
end
H = sqrt(R_tunnel^2 - X.^2) + H_INIT;
    
Z_arr = linspace(0,L_total,2*N);
Z = Z_arr(1)*ones(1,N);
for i = 1:2*N-1
    Z = [Z Z_arr(i+1)*ones(1,N)];
end

scatter3(X,Z,H,'.');
xlim = ([-0.1 0.1]);
ylim = ([0 0.2]);


%% Finding amount of sand + tungst
L_box = 12*0.0254; % Length of box or tunnel section
H_box = 11*0.0254; % Total box height
W_box = 8*0.0254; % Width of box
vol_tunnel = pi*R_tunnel^2*L_box; % Volume of full 3-ring segment

fill_height = 10*0.0254; % Height to which material filled
prct_sand = 0.90; prct_tungst = 1 - prct_sand;
Vol_soil = L_box*fill_height*W_box - vol_tunnel; 
Vol_Sand = Vol_soil*prct_sand;
Vol_Tungst = Vol_soil-Vol_Sand;
H_Tungst = Vol_Tungst/(L_box*W_box);

m_Sand = Vol_Sand*rho_sand;
m_Tungst = Vol_Tungst*rho_tungst;


%% Calculating possible pressure applied by sand

H_mid_to_top = 7*0.0254; % Height of middle of tunnel to top of box 
H_top_to_top = 4*0.0254; % Height of top of tunnel to top of box

rho_total = prct_sand*rho_sand + prct_tungst*rho_tungst; % kg/m^3

vol_box_min = H_top_to_top * W_box * L_box; % The box directly above the tunnel
vol_box_mid = H_mid_to_top * W_box * L_box; % The box from mid of tunnel to top of box
vol_tunnel_half = vol_tunnel/2;
vol_sand_min = vol_box_min;
vol_sand_mid = vol_box_mid - vol_tunnel_half;

mass_min = vol_sand_min * rho_total;
mass_mid = vol_sand_mid * rho_total;

%% Testing the model with centrifuge speeds
RPM_arr = [0:60]; % Testing for RPM speeds ranging from 0 to 60
R_cent = 1.5; % Radius of centrifuge arm in meters
V_arr = (2*pi*R_cent/60)*RPM_arr;

F_arr = mass_min * V_arr.^2/R_cent;
Area = 2*R_tunnel*L_total;
P_array = F_arr/Area;
% P_array = P_array*2;

acc_arr = F_arr/mass_min;

%% Hoop strain

% Test Tunnel segment dimensions
E_tarp = 1.07*10^9; % Youngs Modulus of tarp, Pa
L_total = 12*2.54/100; % Length of whole segment in meters
L_rings = L_total/2; % Separation center-to-center between two rings
R_tunnel = 3*2.54/100; % Radius of a ring, calculated from rod length [m]
Circ_tunnel = 2*pi*R_tunnel; % Circumference of tunnel in meters

t_tarp = 1e-3; % Thickness of tarp in meters

% Pressure data analysis
P_array_2 = -(0:0.1:4.5); % Array of pressures to test [psi]
P_array_2 = P_array_2 * 6894.76; % Converting to Pa

% Calculating stresses and strain
sig_1 = P_array * R_tunnel / t_tarp; % Across surface of tunnel
sig_2 = P_array * R_tunnel / (2 * t_tarp); % Down Length of tunnel
strain = sig_2/E_tarp; 
def = strain * L_rings; % How much the tunnel would stretch axially 


%% Using straight up pressure on a plane
equiv_width = 1*0.0254;
equiv_area = equiv_width * 1;
I = 1/12 * equiv_width * (t_tarp)^3;
w = P_array * equiv_area; % W in N/m
def_new = 4*w*L_rings^4/(384*E_tarp*I);

%% What about just E = stress/strain

Strain = P_array/E_tarp;
def_3 = Strain * L_rings;

%% centroid
m_box = 22.055; %lbs
m_sand = 51.105-23.765; %lbs
m_tungst = 64-42; %lbs
M_TOTAL = m_box+m_sand+m_tungst;

y_box = H_box/2;
y_sand = H_mid_to_top/2 + H_mid_to_top;
y_tungst = H_box-(1.5*0.0254);

centroid = (m_box*y_box + m_sand*y_sand + m_tungst*y_tungst)/(M_TOTAL);
centroid_in = centroid/0.0254;

%% moments

L_cent = (53.5-centroid_in);
M_CTRW = L_cent*M_TOTAL/47.5;

L_CTRW = [40:0.1:53.5];
M_ARRAY = L_cent*M_TOTAL./L_CTRW;

figure;
plot(53.5-L_CTRW, M_ARRAY);

 
%% Final calcs;
Max_G = 18;
Max_accel = Max_G*9.81; %m/s^2
Frac_sand = H_mid_to_top/(H_box-2.5*0.0254);
M_TOP = (m_sand + m_tungst);
M_TOP_k = convmass(M_TOP,'lbm','kg');

F_FINAL = M_TOP_k*Max_accel;
P_FINAL = F_FINAL/(L_box*W_box);
FOS = P_FINAL/35000;

rho_powd = convmass(m_tungst,'lbm','kg')/(L_box*W_box*0.0254);
Vol_sand = (L_box*W_box*(11-2.5)*0.0254) - (pi*R_tunnel^2*L_box);
rho_sand_act = convmass(m_tungst,'lbm','kg') / Vol_sand;
