% ASEN 2004 Milestone 2 (Lab 1)
% Glider Design:
%Date of Submission: February 28th,2021
% Bennett Grow, Graham Kersey, Zak Reichenbach, 
% Alfredo Bautista, Ketan Kamat, Jackson DePenning

clc 
clear
close all
g = 9.81;

% Control what plots to display here:
plot_lift_curves = true;
plot_drag_polar = true;

%% Dimensions of Glider:
% All center of gravity measurements taken from the leading edge of wing
% Measurements in front of wing have negative values

x_end = 0.563; % Distance from wing leading edge to end of airplane
dihedral = 2.2; % [deg]

% Wing
AR = 3;
wingspan = .62;
c = wingspan/AR; % chord length
Sref = c * wingspan;
W_wing = 2 * Sref * (0.295*9.81);
xCG_wing = c / 2;
CG_wing = xCG_wing * W_wing;

% Vertical Stabilizer
ver_c = 0.07;
ver_span = 0.1;
Sv = ver_c * ver_span; % vertical tail planform area
AR_ver = ver_span^2 / Sv;
xCG_ver = x_end - (2 * ver_c);
W_ver = Sv * (0.295*9.81);
CG_ver = xCG_ver * W_ver;

% Horizontal Stabilizer
hor_c = 0.07;
hor_span =  0.32;
Sh = hor_c * hor_span; % horizontal tail planform area
AR_hor = hor_span^2 / Sh;
xCG_hor = x_end - (2 * hor_c);
W_hor = Sh * (0.295*9.81);
CG_hor = xCG_hor * W_hor;

% Fuselage
fuselage_area = 0.12;
xCG_fuselage = 0.201; % Likely not accurate
W_fuselage = fuselage_area * (0.295*9.81);
CG_fuselage = xCG_fuselage * W_fuselage;

% Payload
m_garmin = 0.16  *g;
xCG_payload = 0.07;
CG_payload = m_garmin * xCG_payload;

% Ballast
m_ballast = 0.05  *g;
xCG_ballast = -0.05;
CG_ballast = xCG_ballast * m_ballast;

% Whole aircraft
Swet = fuselage_area + 2*(Sref + Sh + Sv);
area_total = Sref + Sh + Sv + fuselage_area;

% Weight and Mass
W = W_wing + W_ver + W_hor + W_fuselage + m_garmin + m_ballast;
m = W/g;

cost = (area_total)*(0.295) * 1000;

% Center of gravity location
xCG = (CG_wing + CG_ver + CG_hor + CG_fuselage + CG_payload + CG_ballast) / W;

wingload = W/Sref;

%% Load Data:
load app2data
alpha = app2data(:,1);
Cl = app2data(:,2);
Cd = app2data(:,3);
Cm4 = app2data(:,4);
Clmax = max(Cl); %Maximum Cl

%% Creating a 3D Lift Curve:
%Creating a Line Fit
alpharange = 9:23;
[p,~] = polyfit(alpha(alpharange), Cl(alpharange),1);
a0 = p(1);
e = 0.8;

% Wing
a = a0/(1+((57.3*a0)/(pi*e*AR)));
CL = a * (alpha(alpharange) - alpha(16));

% Horizontal Stab.
a_hor = a0/(1+((57.3*a0)/(pi*e*AR_hor)));
CL_hor = a_hor * (alpha(alpharange) - alpha(16));

% Vertical Stab.
a_ver = a0/(1+((57.3*a0)/(pi*e*AR_ver)));
CL_ver = a_ver * (alpha(alpharange) - alpha(16));


%% Constants
rho = 1.0581;
Cfe = 0.003; %Skin Friction
e0 = 1.78*(1-0.045*AR^0.68)-0.64; %Oswald Efficiency 
k = 1/(pi*e*AR);
k0 = 1/(pi*e0*AR);
h = 7; %Launch Height (m)

%% Drag
CD0 = Cfe*(Swet./Sref); %Parasite Drag Coefficient
CD_wing = Cd(alpharange) + (k * CL.^2);
CD_whole = CD0 + (k0* CL.^2);

%% Glide Range
CL_r = sqrt(CD0/k); %Coefficient of Lift for Glide Range
glide_range = (sqrt(CD0/k)./(2*CD0))*h; %Glide Range (m)

%% Max Endurance
CD_e = 4.*CD0; %Coefficient of Drag for Endurance
Ve = sqrt(wingload./(0.5.*rho.*sqrt(3.*CD0.*pi.*e0.*AR))); % Endurance Velocity (m/s)
qe = 0.5.*rho.*Ve.^2; %Dynamic Pressure for Endurance (Pa)
E = wingload.*h./(CD_e.*qe.*Ve); % Endurance (s)
CL_e = sqrt(3.*CD0*pi*e0*AR); %Coefficient of Lift for Glide Endurance

%% Velocity Required for Max E and R
Vr = sqrt(wingload./(sqrt(CD0.*pi.*e0.*AR).*0.5.*rho)); % Range Velocity (m/s)
Vstall = sqrt(2*W./(Clmax.*rho.*Sref)); %Stall Velocity (m/s)

%% Stability
% Horizontal tail volume coef
VH = Sh*((xCG_hor + 0.25 * hor_c) - xCG)/(Sref * c);

% % CL at max range or endurace, use whole aircraft CL
CLw = CL_r;

% Calculate angle of attack where CL of wing = CLr
[~, AOAindex] = min(abs(CL-CL_r));
AOA = alpha(alpharange(1) + AOAindex - 1);

% Coef of moment about aero center of wing, approx w/ 2D airfoil value
Cmacw = Cm4(AOA);

% CL for horizontal tail
it = 0; %SOLVE FOR THIS
At = AOA -it;
CLht = a_hor*AOA + a_hor*it; %CL_hor(AOAindex)

% Actual pitching moment
Cmcg = Cmacw + CLw*(0.25 * c)-VH*CLht;

% Ideal CLht (set Cmcg=0)
CLht_ideal = (Cmacw + CLw * 0.25 * c)/VH;


% Tail angle of attack
%    Assumes downwash is 0
alpha_tail = AOA - it;

% Vertical tail volume coef
VV = (Sv * (xCG_ver - 0.25 * ver_c))/(Sref * wingspan);

% Spiral parameter B
B = ((xCG_ver - 0.25 * ver_c) * dihedral) / (wingspan * CLw);

% Neutral Point
Hnp = .25 + VH*(a_hor/a);
Xnp = (Hnp * c);

% Static Margin
SM = Hnp - (xCG/c);

%% Trim Diagrams:

% Maximum Range:
alphastall = alpha(28);
Cmacw_trim = Cm4(AOA);
AOAtrim_r = alpha(16:36);
CLht_trim_r = a_hor*AOAtrim_r +a_hor*it;
Cmcg_trim_r = Cmacw_trim + CL_r*(xCG/c-0.25)-VH*CLht_trim_r;

%Plotting Trim Diagram:
plot(AOAtrim_r,Cmcg_trim_r)
title(" Trim Graph for Maximum Range");
xline(alphastall,'--','Color','m');
xlabel("Angle of Attack (deg)");
ylabel("Cmcg");
yline(0);
grid on
legend("Trim Data","alpha-stall");

%Set Endurance Range:
alphastall = alpha(28);
Cmacw_trim = Cm4(AOA);
AOAtrim_e =alpha(16:36);
CLht_trim_e = a_hor*AOAtrim_e +a_hor*it;
Cmcg_trim_e = Cmacw_trim + CL_e*(xCG/c-0.25)-VH*CLht_trim_e;

%Plotting Trim Diagram:
figure
plot(AOAtrim_e,Cmcg_trim_e)
title("Trim Diagram for Max Endurance");
xline(alphastall,'--');
xlabel("Angle of Attack (deg)");
ylabel("Cmcg");
yline(0);
grid on
xline(alphastall,'--','Color','m');
xlabel("Angle of Attack (deg)");
ylabel("Cmcg");
yline(0);
grid on
legend("Trim Data","alpha-stall");
%% Print calculations
fprintf("Max Glide Range:                %.3f m  \n", glide_range);
fprintf("Range Velocity:                 %.3f m/s  \n", Vr);
fprintf("Max Glide Endurance:            %.3f sec  \n", E);
fprintf("Endurance Velocity:             %.3f m/s  \n", Ve);
fprintf("Cost:                           $%.2f  \n", cost);
fprintf("CG Percentage:                  %.f%% of chord  \n", xCG/c * 100);
fprintf("Long. Static Stability V_H:     %.3f  \n", VH);
fprintf("Lateral Static Stability V_V:   %.3f  \n", VV);
fprintf("Spiral Parameter B:             %.2f with %.2f deg. dihedral \n", B, dihedral);
fprintf("Wing Loading:                   %.2f  \n", wingload);
fprintf("Coef. of Moment about CG:       %.2f   \n", Cmcg);
fprintf("Neutral Point Distance:         %.2f  \n",Xnp);
fprintf("Static Margin:                  %.2f   \n",SM);

%% Plotting

% Lift Curves
if plot_lift_curves == true
    figure
    hold on
    grid on
    
    plot(alpha(6:31), Cl(6:31), 'LineWidth', 2)
    plot(alpha(alpharange), CL, 'LineWidth', 2)
    plot(alpha(alpharange), CL_ver, ':', 'LineWidth', 2)
    plot(alpha(alpharange), CL_hor, ':', 'LineWidth', 2)

    yline(0);
    xline(0);
    xlabel('Angle of Attack (Deg)')
    ylabel('Coefficient of Lift')
    title('Wing 2D and 3D Lift Curves')
    legend('2D','3D Wing', '3D Vert. Stab.', '3D Hor. Stab.', 'Location','SouthEast')
    xlim([-10 10])
    ylim([-0.8 0.8])
    hold off
end

% Drag Polar
if plot_drag_polar == true
    figure
    hold on
    grid on
    
    plot(CL, Cd(alpharange), 'LineWidth', 2);
    plot(CL, CD_wing, 'LineWidth', 2)
    plot(CL, CD_whole, 'LineWidth', 2)
    
    yline(0);
    xline(0);
    legend('2D Airfoil Drag Polar', '3D Wing Drag Polar','Whole Aircraft Drag Polar','Location','North')
    xlabel('Coefficient of Lift')
    ylabel('Coefficient of Drag')
    title('Drag Polar')
    hold off
end









