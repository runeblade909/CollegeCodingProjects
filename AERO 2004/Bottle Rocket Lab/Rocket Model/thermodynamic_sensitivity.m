%{
thermodynamic_sensitivity.m
ASEN2004 Spring 2021 Bottle Rocket Lab
Author: Bennett Grow

Calculates the thrust of a water bottle rocket based on the thermodynamics 
of water expulsion and air expansion then models the motion of flight.

Essentially the same model as used in Project 2 of ASEN2012, but uses
vector instead of trig and models in 3D instead of 2D. 

This version performs sensitivity analyses on 5 parameters:
    - Coefficent of drag
    - Mass of water propellant
    - Density of the water propellant
    - Temperature of the water propellant
    - Launch pad angle

Requires function 'water_density.m'

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

%% Loop Setup
% Housekeeping
clc;
clear;
close all;

% Plot graphs?
enable_plots = true;
% Number of variables to perform sensitivity analysis on
num_sensitivity = 5;
% Number of trials per sensitivity variable
num_trials = 20;
% Total number of simualtions to run
simulations = num_sensitivity * num_trials;

% Preallocate for to carry simulation data outside of the k for loop
data = cell(num_sensitivity, num_trials);           % States
dataxy = cell(num_sensitivity, num_trials);         % Normalized xy vectors
sens = zeros(num_sensitivity, num_trials);          % Sensitivity values
coorx = zeros(num_sensitivity, num_trials);         % End x coordinates
coory = zeros(num_sensitivity, num_trials);         % End y coordinates
coorxy = zeros(num_sensitivity, num_trials);        % End distance
maxz = zeros(num_sensitivity, num_trials);          % Max trajectory heights

for i = 1:num_sensitivity
    for k = 1:num_trials

    s = struct; % Struct to hold all info to ease passing to ODE function

    %% Launch Day & Sensitivity Analysis Variables
    s.stand_ang         =       45;          %[deg] Test stand angle
    s.P_gage            =       40;          %[Psi] Initial bottle gage pressure  
    s.C_D               =       0.4;         % Coef. of drag (0.3-0.5)
    air_temp            =       63;          %[deg F]
    prop_mass           =       600;         %[g]
    dry_mass            =       160;         %[g]
    wind_speed          =       8;           %[mph]
    wind_direction      =       225;         %[deg clockwise from N]
    launch_heading      =       40;          %[deg clockwise from N]
    s.rho_w             =       1000;        %[kg/m^3] Density of water

    %% Sensitivity Ranges
    % Each range will be the +/- added to the sensitivity variables
    coef_drag_range     =       0.1;         %
    prop_mass_range     =       100;         %[g]
    prop_density_range  =       10;         %[kg/m^3]
    stand_ang_range     =       10;          %[deg]

    % Create num_trials of linearly spaced values for each sensitivity range 
    delta_coef_drag     = linspace(-coef_drag_range, coef_drag_range, num_trials);
    delta_prop_mass     = linspace(-prop_mass_range, prop_mass_range, num_trials);
    delta_prop_density  = linspace(-prop_density_range, prop_density_range, num_trials);
    delta_prop_temp     = linspace(35, 90, num_trials);         %[deg F]
    delta_stand_ang     = linspace(-stand_ang_range, stand_ang_range, num_trials);

    % Choose which sensitivity analysis to perform
    % First num_trials will be coef_drag, second set of num_trials will be prop mass, etc.
    if (i == 1) % Coef of drag
        s.C_D = s.C_D + delta_coef_drag(k);
        sens(i, k) = s.C_D;
    elseif(i == 2) % prop mass
        prop_mass = prop_mass + delta_prop_mass(k);
        sens(i, k) = prop_mass;
    elseif(i == 3) % prop density
        s.rho_w = s.rho_w + delta_prop_density(k);
        sens(i, k) = s.rho_w;
    elseif(i == 4) % prop temp
        s.rho_w = water_density(delta_prop_temp(k), 'F');
        sens(i, k) = delta_prop_temp(k);
    elseif(i == 5) % stand angle
        s.stand_ang = s.stand_ang + delta_stand_ang(k);
        sens(i, k) = s.stand_ang;
    end


    %% Universal Constants
    s.g           =       -9.81;                                                %[m/s^2]
    s.gamma       =       1.4;                                                  % Ratio of specific heats for air
    s.R           =       287;                                                  %[J/kg K] Uni. Gas Const.
    s.rho_amb     =       0.961;                                                %[kg] Ambient air density
    s.P_amb       =       83426.56;                                             %[Pa] (12.1 psi) Ambient air pressure

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
    z = 0.25;
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

    for j = 1:length(state)
        [~, m_r(j), heading(j, :), F(j, 1:3)] = odefunc(t, state(j, :), s);
        F(j, 4) = sqrt(F(j, 1)^2 + F(j, 2)^2 + F(j, 3)^2);
        xy(j) = sqrt(x(j)^2 + y(j)^2);
    end

    %% Combine this run's info into pre-loop allocated vectors
    data{i, k} = state;
    dataxy{i, k} = xy;
    coorx(i, k) = x(end);
    coory(i, k) = y(end);
    coorxy(i, k) = xy(end);
    maxz(i, k) = max(z);

    end
end

%% Calculate Max Values and Print
[C_D_max, C_D_i] = max(coorxy(1,:));
[prop_mass_max, prop_mass_i] = max(coorxy(2,:));
[prop_density_max, prop_density_i] = max(coorxy(3,:));
[prop_temp_max, prop_temp_i] = max(coorxy(4,:));
[stand_ang_max, stand_ang_i] = max(coorxy(5,:));

fprintf("Coef. of Drag:     %.2f,     Max Range: %.2f \n", sens(1, C_D_i), C_D_max);
fprintf("Prop. Mass:        %.2f,   Max Range: %.2f \n", sens(2, prop_mass_i), prop_mass_max);
fprintf("Prop. Density:     %.2f,   Max Range: %.2f \n", sens(3, prop_density_i), prop_density_max);
fprintf("Prop. Temp:        %.2f,    Max Range: %.2f \n", sens(4, prop_temp_i), prop_temp_max);
fprintf("Stand Angle:       %.2f,    Max Range: %.2f \n", sens(5, stand_ang_i), stand_ang_max);

%% Plotting
if enable_plots == true
    
    %--- Variation Plotting ---%
    
    % Coef. of Drag
    figure(1);
    plot(coorxy(1,:), sens(1,:), 'bo');
    title('Coefficent of Drag Sensitivity Analysis');
    ylabel('Distance [m]');
    xlabel('Coefficent of Drag');
    
    % Prop mass
    figure(2);
    plot(coorxy(2,:), sens(2,:), 'bo');
    title('Propellant Mass Sensitivity Analysis');
    ylabel('Distance [m]');
    xlabel('Propellant Mass [g]');
    
    % Prop density
    figure(3);
    plot(coorxy(3,:), sens(3,:), 'bo');
    title('Propellant Density Sensitivity Analysis');
    ylabel('Distance [m]');
    xlabel('Propellant Density [kg/m^3]');
    
    % Prop temp
    figure(4);
    plot(coorxy(4,:), sens(4,:), 'bo');
    title('Propellant Temperature Sensitivity Analysis');
    ylabel('Distance [m]');
    xlabel('Propellant Temperature [F]');
    
    % Stand ang
    figure(5);
    plot(coorxy(5,:), sens(5,:), 'bo');
    title('Launch Angle Sensitivity Analysis');
    ylabel('Distance [m]');
    xlabel('Launch Angle [deg]');
    
    %--- Trajectory Plotting ---%
    %{
    % Convert sens to string matrix for legends
    sens_string = arrayfun(@num2str, sens, 'un', 0);

    % Coef. of drag
    figure(1);
    hold on;
    axis equal;
    for i = 1:num_trials
        plot(dataxy{1,i}, data{1,i}(:,3));
    end
    xlim([0 coorxy(1, C_D_i)]);
    ylim([0 30]);
    title("Coefficent of Drag Sensitivity Analysis");
    xlabel('Distance (xy) [m]'); 
    ylabel('Height (z) [m]');
    lgd = legend(sens_string(1,:));
    lgd.Orientation = 'horizontal';
    lgd.NumColumns = 4;
    lgd.Location = 'southoutside';
    grid on;
    hold off;

    % Prop Mass
    figure(2);
    hold on;
    axis equal;
    for i = 1:num_trials
        plot(dataxy{2,i}, data{2,i}(:,3));
    end
    xlim([0 coorxy(2, prop_mass_i)]);
    ylim([0 30]);
    title("Propellant Mass Sensitivity Analysis");
    xlabel('Distance (xy) [m]'); 
    ylabel('Height (z) [m]');
    lgd = legend(sens_string(2,:));
    lgd.Orientation = 'horizontal';
    lgd.NumColumns = 4;
    lgd.Location = 'southoutside';
    grid on;
    hold off;

    % Prop Density
    figure(3);
    hold on;
    axis equal;
    for i = 1:num_trials
        plot(dataxy{3,i}, data{3,i}(:,3));
    end
    xlim([0 coorxy(3, prop_density_i)]);
    ylim([0 30]);
    title("Propellant Density Sensitivity Analysis");
    xlabel('Distance (xy) [m]'); 
    ylabel('Height (z) [m]');
    lgd = legend(sens_string(3,:));grid on;
    lgd.Orientation = 'horizontal';
    lgd.NumColumns = 4;
    lgd.Location = 'southoutside';
    hold off;

    % Prop Temp
    figure(4);
    hold on;
    axis equal;
    for i = 1:num_trials
        plot(dataxy{4,i}, data{4,i}(:,3));
    end
    xlim([0 coorxy(4, prop_temp_i)]);
    ylim([0 30]);
    title("Propellant Temperature Sensitivity Analysis");
    xlabel('Distance (xy) [m]'); 
    ylabel('Height (z) [m]');
    lgd = legend(sens_string(4,:));grid on;
    lgd.Orientation = 'horizontal';
    lgd.NumColumns = 4;
    lgd.Location = 'southoutside';
    hold off;

    % Stand Angle
    figure(5);
    hold on;
    axis equal;
    for i = 1:num_trials
        plot(dataxy{5,i}, data{5,i}(:,3));
    end
    xlim([0 coorxy(5, stand_ang_i)]);
    ylim([0 30]);
    title("Launch Stand Angle Sensitivity Analysis");
    xlabel('Distance (xy) [m]'); 
    ylabel('Height (z) [m]');
    lgd = legend(sens_string(5,:));grid on;
    lgd.Orientation = 'horizontal';
    lgd.NumColumns = 4;
    lgd.Location = 'southoutside';
    hold off;
    %}
end

%% ODE Integration Function
function [state, m_r, heading, F] = odefunc(t, state0, s)

    % Extract integrated vars from state vector
    x       = state0(1);
    y       = state0(2);
    z       = state0(3);
    vx      = state0(4);
    vy      = state0(5);
    vz      = state0(6);
    m_air   = state0(7);
    vol_air = state0(8);

    s.rho_w
    
    
    % Total rocket mass
    m_r = s.m_B + s.rho_w*(s.vol_B - vol_air) + m_air
    
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