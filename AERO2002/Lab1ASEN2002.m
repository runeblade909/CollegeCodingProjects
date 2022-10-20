 %% High Altitude Balloon Calculator
% Written by Bennett Grow, Amanda Marlow, & Zak Reichenbach 9/2020

% Finds the mass of gasses and the radius of the balloon required to
% fly a 500 kg payload at 25km

% Assumptions:
%   Balloon is spherical
%   No transient effects during ascent
%   No stresses from payload or wind loading
%   Payload mass includes all required support structure

clear
clc

%% US Standard Atmosphere Air Properties at 25km from Eng. Toolbox:
altitude_vec = 20000 : 10 : 27000;
[std_T, std_SOS, std_P, std_Rho] = atmoscoesa(altitude_vec);
[std_T_26, std_SOS_26, std_P_26, std_rho_26] = atmoscoesa(26000);
[std_T_25, std_SOS_25, std_P_25, std_rho_25] = atmoscoesa(25000);
[std_T_24, std_SOS_24, std_P_24, std_rho_24] = atmoscoesa(24000);

%% Day and night temp due to thermal radiation
alpha_sb = 0.6;
alpha_eb = 0.8;
epsilon_b = 0.8;
alpha_SB = 5.670 * 10^(-8); %[J/K^4*m^2*s]
q_sun = 1353; %[W/m^2]
q_earth = 237; %[W/m^2]

T_night = nthroot((alpha_eb * q_earth)/(4 * epsilon_b * alpha_SB), 4);
T_day = nthroot((alpha_eb * q_earth + alpha_sb * q_sun)/(4 * epsilon_b * alpha_SB), 4);

%% Constants:
FOS = 1.75;
P_g = 10; % [Pa] Gauge pressure
R = 8.31446; % [M^3*Pa/(mol*K)] Universal Gas Constant
payload_mass = 500; % [kg]

%% Gas properties
Mh2 = 2.016/1000; % [kg/mol] molar mass of hydrogen gas
density_h2_day_25 = (std_P_25 + P_g) / ((R/Mh2) * T_day); % [g/L] or [kg/m^3] Ideal gas law
density_h2_night_25 = (std_P_25 + P_g) / ((R/Mh2) * T_night);
density_h2_day_26 = (std_P_26 + P_g) / ((R/Mh2) * T_day);

%% Medium density polyethylene (MDPE), film grade properties:
ys = 19.9 * 10^6; % [Pa] Yield strength
density_bal = 937; % [kg/m^3] density
% http://www.matweb.com/search/DataSheet.aspx?MatGUID=d9e80589aaf144a3a8eaad08af77f181&ckck=1

%% Solve for radius R in boyancy eqn
syms r;
a = (4/3)*pi;
r_ext =  r + (FOS * P_g * r)/(2*ys);
radius_day_eqn = std_rho_25 * a * r_ext^3 == payload_mass + density_bal * a * (r_ext^3 - r^3) + density_h2_day_25 * a * r^3;
radius_day_all = vpasolve(radius_day_eqn, r);

radius_day = double(radius_day_all(1));
radius_day_ext =  radius_day + (FOS * P_g * radius_day)/(2*ys);

%% Calculate thickness, volume, and weight then print
thickness = (FOS * P_g * radius_day)/(2 * ys) *1000000; % [microns] Material thickness eqn from lab videos
volume_day = double(a * radius_day_ext^3);
mass_h2_day = volume_day * density_h2_day_25;
sys_mass_day = double( payload_mass + density_bal * a * (radius_day_ext^3 - radius_day^3) + mass_h2_day);
sys_density_day = double(sys_mass_day/volume_day);

% Does the mass need to be vpa solved with ballast removed???
volume_night = mass_h2_day/density_h2_night_25 ;
sys_density_night = mass_h2_day / volume_night;
sys_mass_night_ballast_dropped = std_rho_24 * volume_night;

syms delta_h2;
vented_h2_eqn = (mass_h2_day - delta_h2) / (sys_mass_night_ballast_dropped - delta_h2) == density_h2_day_26 / std_rho_26;
vented_h2 = vpasolve(vented_h2_eqn, delta_h2);

ballast = sys_mass_day - sys_mass_night_ballast_dropped;

volume_day_2 = (mass_h2_day - vented_h2) / density_h2_day_26;

sys_density_day_2 = (sys_mass_night_ballast_dropped - vented_h2) / volume_day_2;




%%

fprintf("1: Day,   25km: %.3fkg total mass,  %.3fkg of H2, %.3fm radius, %.3f micron thick MDPE",  sys_mass_day, mass_h2_day, radius_day, thickness);
fprintf("\n2: Night, 25km: %.3fkg total mass,  %.3fkg of H2", sys_mass_night_ballast_dropped, mass_h2_day);
fprintf("\n3: Day 2, 25km: %.3fkg total mass,  %.3fkg of H2", (sys_mass_night_ballast_dropped - vented_h2), (mass_h2_day - vented_h2));

fprintf("\nWe dropped %.3fkg of ballast mass and vented out %.3fkg of hydrogen to maintain the 25km +- 1km range.\n", ballast, vented_h2);


%%

















