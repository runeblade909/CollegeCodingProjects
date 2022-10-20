%{
thermodynamic_montecarlo.m
ASEN2004 Spring 2021 Bottle Rocket Lab
Author: Bennett Grow

Calculates the thrust of a water bottle rocket based on the thermodynamics 
of water expulsion and air expansion then models the motion of flight.

Essentially the same model as used in Project 2 of ASEN2012, but uses
vector instead of trig and models in 3D instead of 2D. 

This

%}

%% Assumptions
%{
- Ignores static test stand data
- Ignores water and air viscosity
- Ignores change in bottle shape and volume due to internal pressure
- Ideal nozzle
- Wind is a constant vector in x and y

x is direction of launch
z is upwards
y is RHR

%}

clc;
clear;
close all;

simulations = 100;

% Preallocate
data = cell(3, simulations);
coorx = zeros(1, simulations);
coory = zeros(1, simulations);
coorxy = zeros(1, simulations);
maxz = zeros(1, simulations);

for k = 1:simulations

s = struct; % Struct to hold all info to ease passing to ODE function

%% Launch Day Variables
s.stand_ang         =       41;          %[deg] Test stand angle
s.P_gage            =       40;          %[Psi] Initial bottle gage pressure  
s.C_D               =       0.3;         % Coef. of drag (0.3-0.5)
air_temp            =       63;          %[deg F]
prop_mass           =       600;         %[g]
dry_mass            =       160;         %[g]
wind_speed          =       8;           %[mph]
wind_direction      =       225;         %[deg clockwise from N]
launch_heading      =       40;          %[deg clockwise from N]

%% Uncertainties 
linvar = 100;
sigma_stand_ang = linspace(-1, 1, linvar);
sigma_water = linspace(-5, 5, linvar);
sigma_C_D = linspace(-0.1, 0.1, linvar);
sigma_wind_speed = linspace(-2.5, 2.5, linvar);
sigma_wind_direction = linspace(-22.5, 22.5, linvar);
sigma_heading = linspace(-3, 3, linvar);

%% Random Uncertainty Selection
pick = randi([1 linvar], 1, 1);
delta_stand_ang = sigma_stand_ang(pick);
pick = randi([1 linvar], 1, 1);
delta_water = sigma_water(pick);
pick = randi([1 linvar], 1, 1);
delta_C_D = sigma_C_D(pick);
pick = randi([1 linvar], 1, 1);
delta_wind_speed = sigma_wind_speed(pick);
pick = randi([1 linvar], 1, 1);
delta_wind_direction = sigma_wind_direction(pick);
pick = randi([1 linvar], 1, 1);
delta_heading = sigma_heading(pick);

%% Apply Randomly Picked Uncertainty
s.stand_ang         =       s.stand_ang + delta_stand_ang;       
s.C_D               =       s.C_D + delta_C_D;                    
prop_mass           =       prop_mass + delta_water;                
wind_speed          =       wind_speed + delta_wind_speed;
wind_direction      =       wind_direction + delta_wind_direction;  
launch_heading      =       launch_heading + delta_heading;       

%% Universal Constants
s.g           =       -9.81;                                                %[m/s^2]
s.gamma       =       1.4;                                                  % Ratio of specific heats for air
s.R           =       287;                                                  %[J/kg K] Uni. Gas Const.
s.rho_amb     =       0.961;                                                %[kg] Ambient air density
s.P_amb       =       83426.56;                                             %[Pa] (12.1 psi) Ambient air pressure
s.rho_w       =       1000;                                                 %[kg/m^3] Density of water

%% Rocket Specifications & Conversions
s.dia_B       =       0.105;                                                %[m] Bottle diameter
s.dia_t       =       0.021;                                                %[m] Throat diameter
s.vol_B       =       0.002;                                                %[m^3] volume of bottle
s.m_B         =       dry_mass / 1000;                                      %[g to kg] Mass of empty bottle

s.A_B         =       0.25*pi*(s.dia_B^2);                                  %[m^2] Crosssectional area of front of bottle
s.A_t         =       0.25*pi*(s.dia_t^2);                                  %[m^2] Throat area

s.C_d         =       0.8;                                                  % Discharge coef.
s.stand_l     =       0.61;                                                 %[m] Test stand length
s.stand_h     =       0.25;                                                 %[m] Rocket starting height on test stand

w_mass_i      =       prop_mass;                                            %[g]
s.vol_w_i     =       w_mass_i * 10^-6;                                     %[g to m^3] Initial volume of water
s.T_air_i     =       (air_temp - 32) * (5/9)  +273.15;                     %[F to K] Initial temp of air

s.P_B_i       =       344738 + s.P_amb;                                     %[Pa] Initial bottle gage pressure (50 psi)
s.vol_air_i   =       s.vol_B - s.vol_w_i;                                  %[m^3]
s.m_air_i     =       ((s.P_B_i) * s.vol_air_i)/(s.R * s.T_air_i);          %[kg]
s.m_r_i       =        s.m_B + s.rho_w * s.vol_w_i + s.m_air_i;             %[kg] initial total mass of rocket
s.P_B_i       =       s.P_gage*6895 + s.P_amb;                              %[Pa] Initial bottle total pressure


%% Motion Information

% Launch Heading Vector
%s.heading_i = [sind(launch_heading)*cosd(s.stand_ang) sind(launch_heading)*sind(s.stand_ang) cosd(s.stand_ang)];
s.heading_i = [cosd(s.stand_ang) 0 sind(s.stand_ang)];

% Wind Vector
wind_speed = wind_speed * 0.44704; %[mph to m/s]
s.v_wind = wind_speed * [sind(wind_direction - launch_heading) cosd(wind_direction - launch_heading) 0];

x = 0;
y = 0;
z = 0;
vx = 0;
vy = 0;
vz = 0;

state0 = [x y z vx vy vz s.m_air_i s.vol_air_i]; % Integrating variables

%% ODE Call
tspan = [0 5];
opt = odeset('Events', @groundstrike); % Stop when the bottle hits the ground
[t, state] = ode45(@(t, state) odefunc(t, state, s), tspan, state0, opt);

x       = state(:, 1);
y       = state(:, 2);
z       = state(:, 3);
vx      = state(:, 4);
vy      = state(:, 5);
vz      = state(:, 6);
m_air   = state(:, 7);
vol_air = state(:, 8);

m_r = zeros(length(state),1);
heading = zeros(length(state),3);
F = zeros(length(state),4);
xy = zeros(length(state),1);

for i = 1:length(state)
    [~, m_r(i), heading(i, :), F(i, 1:3)] = odefunc(t, state(i, :), s);
    F(i, 4) = sqrt(F(i, 1)^2 + F(i, 2)^2 + F(i, 3)^2);
    xy(i) = sqrt(x(i)^2 + y(i)^2);
end

%% Combine this run's info
data{1, k} = state;
data{2, k} = F;
data{3, k} = xy;

coorx(k) = x(end);
coory(k) = y(end);
coorxy(k) = xy(end);
maxz = max(z);

end

drift_angle = atand(mean(coory)/mean(coorx));
fprintf("Averages: Max Height: %.2f m, Distance: %.2f m, Endpoint: (%.2f, %.2f), Drift: %.2f deg \n", mean(maxz), mean(coorxy), mean(coorx), mean(coory), drift_angle);

%% Error Ellipses

figure; plot(coorx,coory,'k.','markersize',6)
axis equal; grid on; xlabel('Downrange (x) [m]'); ylabel('Crossrange (y) [m]'); hold on;
title("Error Ellipses");

% Calculate covariance matrix
P = cov(coorx,coory);
mean_x = mean(coorx);
mean_y = mean(coory);
 
% Calculate the define the error ellipses
n=100;          % Number of points around ellipse
p=0:pi/n:2*pi;  % angles around a circle
 
[eigvec,eigval] = eig(P); % Compute eigen-stuff
xy_vect = [cos(p'),sin(p')] * sqrt(eigval) * eigvec'; % Transformation
x_vect = xy_vect(:,1);
y_vect = xy_vect(:,2);
 
% Plot the error ellipses overlaid on the same figure
plot(1*x_vect+mean_x, 1*y_vect+mean_y, 'b')
plot(2*x_vect+mean_x, 2*y_vect+mean_y, 'g')
plot(3*x_vect+mean_x, 3*y_vect+mean_y, 'r')

hold off

%% Trajectory Plotting

figure;
hold on;
axis equal;
for i = 1:simulations
    plot3(data{1,i}(:,1), data{1,i}(:,2), data{1,i}(:,3), 'Color', '#808080');
end
errorzplane = zeros(length(x_vect), 1);
plot3(1*x_vect+mean_x, 1*y_vect+mean_y, errorzplane, 'b', 'LineWidth', 1);
plot3(2*x_vect+mean_x, 2*y_vect+mean_y, errorzplane, 'g', 'LineWidth', 1);
plot3(3*x_vect+mean_x, 3*y_vect+mean_y, errorzplane, 'r', 'LineWidth', 1);

xlim([0 90]);
ylim([-30 0]);
zlim([0 40]);

xlabel('Downrange (x) [m]'); 
ylabel('Crossrange (y) [m]');
zlabel('Height (z) [m]');
grid on;
view([-1 -1 1]);
hold off;

%% ODE Integration Function
function [state, m_r, heading, F] = odefunc(t, state0, s)
%function [state, m_r] = odefunc(t, state0, s)

    % Extract integrated vars from state vector
    x       = state0(1);
    y       = state0(2);
    z       = state0(3);
    vx      = state0(4);
    vy      = state0(5);
    vz      = state0(6);
    m_air   = state0(7);
    vol_air = state0(8);

    % Total rocket mass
    m_r = s.m_B + s.rho_w*(s.vol_B - vol_air) + m_air;
    
    % Heading Vector
    % If rocket hasn't left test stand
    if (z <= (s.stand_h + s.stand_l*sind(s.stand_ang))) || (x <= s.stand_l*cosd(s.stand_ang)) 
        s.g = 0;
        % Heading
        heading = s.heading_i;
        % Magnitude of velocity
        v = [vx vy vz];
        v_mag = sqrt(vx^2 + vy^2 + vz^2);
        
    else
        s.g = 9.81;
        % Effect of wind
        vx = vx + s.v_wind(1);
        vy = vy + s.v_wind(2);
        vz = vz + s.v_wind(3);
        % Magnitude of velocity
        v = [vx vy vz];
        v_mag = sqrt(vx^2 + vy^2 + vz^2);
        % Heading
        heading = v ./ v_mag;
    end
    
    
    % Phase 1
    % Air pressure exhausts water until there is no more water left in
    % the bottle.
    if vol_air <= s.vol_B
        phase = 1;
        % Change in volume over time
        dvol_air_dt = s.C_d * s.A_t * sqrt((2/s.rho_w)*(s.P_B_i * ((s.vol_air_i ./ vol_air)^s.gamma) - s.P_amb));
        
        % Change in pressure
        P_air = s.P_B_i * (s.vol_air_i ./ vol_air).^s.gamma;
        
        % Force from thrust
        F = 2 * s.C_d * s.A_t * (P_air - s.P_amb);
        
        % Change in air mass
        dm_air_dt = 0; % No air has l eft bottle yet
        
    else % Calculate air p & temp at the end of phase 1
        P_end = s.P_B_i * (s.vol_air_i / s.vol_B)^s.gamma;      % Never changes after phase 1 -  Eqn 15
        T_end = s.T_air_i * (s.vol_air_i / s.vol_B)^s.gamma;    % Never changes after phase 1 -  Eqn 15
        P_air = P_end * (m_air / s.m_air_i)^s.gamma;            % Pressure dependent on air mass after phase 1 - Eqn 14
        

        % Phase 2
        % All water is exhausted. Remaining air is vented until pressure
        % matches ambient pressure. 
        if (P_air > s.P_amb)
            phase = 2;

            % Change in volume over time
            dvol_air_dt = 0;

            % Density and Temp
            roe_air = m_air / s.vol_B;      % Eqn 15
            T = P_air ./ (roe_air .* s.R);  % Eqn 15

            % Critical Pressure
            P_crit = P_air * (2 / (s.gamma + 1))^(s.gamma / (s.gamma - 1));   % Eqn 16

            % 2 cases:
            %   1) p_crit > p_amb  ->  choked flow (M_exit = 1)
            %   2) p_crit < p_amb

            if P_crit > s.P_amb % Choked Flow
                T_exit = (2 / (s.gamma + 1)) * T;                 % Eqn 18
                P_exit = P_crit;                                % Eqn 18
                roe_air_exit = P_exit / (s.R * T_exit);           % Eqn 18
                v_exit = sqrt(s.gamma * s.R * T_exit);              % Eqn 17

            else % p_crit < p_amb
                M_exit = sqrt( ((P_air / s.P_amb)^((s.gamma - 1)/s.gamma) - 1) * (2 / (s.gamma - 1)) );    % Eqn 19
                T_exit = T / (1 + ((s.gamma - 1)/2) * M_exit^2);  % Eqn 20
                roe_air_exit = s.P_amb ./ (s.R .* T_exit);          % Eqn 20
                P_exit = s.P_amb;                                 % Eqn 20
                v_exit = M_exit * sqrt(s.gamma * s.R * T_exit);     % Eqn 21

            end

            % Change in air mass
            dm_air_dt = -s.C_d * roe_air_exit * s.A_t * v_exit;      % Eqn 23

            % Force from thrust
            F = -dm_air_dt * v_exit + (s.P_amb - P_exit) * s.A_t;    % Eqn 22

        % Phase 3
        % Ballistic. No thrust. All water is exhausted and internal pressure
        % matches ambient air pressure. 
        else
            phase = 3;
            F = 0;                                              % Eqn 25
            dm_air_dt = 0;
            dvol_air_dt = 0;

        end
    end

    % *** Calculations that are the same for every phase: *** %
    
    % Drag
    drag = 0.5 * s.rho_amb * v_mag.^2 * s.C_D * s.A_B;
    drag = drag .* heading;
    
    % Force Components
    F = F .* heading;
    
    % Acceleration
    ax = (F(1) - drag(1)) / m_r;
    ay = (F(2) - drag(2)) / m_r;
    az = (F(3) - drag(3) - s.g*m_r) / m_r;
    a = [ax ay az];
    
    % State output to integrate
    state = zeros(8, 1);
    state(1) = vx;
    state(2) = vy;
    state(3) = vz;
    state(4) = ax;
    state(5) = ay;
    state(6) = az;
    state(7) = dm_air_dt;
    state(8) = dvol_air_dt;
    
    output = [phase F vol_air v drag];
    %fprintf("phase %d, m_air: %.2vf, vol_air: %.2f\n", phase, m_air, vol_air);

end

% groundstrike function stops integration when z=0
function [value, isterminal, direction] = groundstrike(t,state)
    value = (state(3) <= 0);
    isterminal = 1;
    direction = 0;
end