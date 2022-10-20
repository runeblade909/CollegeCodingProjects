function [w] = model_4(Theta)
% Purpose: calculations for model 4 - models a unbalanced wheel with a support structure as it
% rolls down a hill. The extra mass is inserted at a certain radius from the
% center, and is modelled as itself with a given radius to account for the angular kinetic energy of the rod mass
% Authors: Jacob Starkel, Krystal Horton, Cate Billings, and Zak Reichenbach
% Date Completed: 3/12/2021

%% Givens 
m_cylinder  = 11.7; % [kg]
m_supports  = 0.7;  % [kg]
m_extra = 3.4; % [kg]
r_cylinder  = .235;  % [m]
k = .203; % [m]
I = m_cylinder*k^2; % [kgm^2]
beta = 5.5; % slope of ramp [degrees]
radius_to_extra = .178; % [m]
g = 9.81; % [m/s^2]
T = [1 1.2 1.4 1.45 1.5]; % torque friction of bearing
m_t = m_cylinder + m_supports; % [kg]

%% Model 4

phi=Theta+(beta*(pi/180));

Num_4 = 2 *(m_t*g*r_cylinder*Theta*sind(beta) + m_extra*g*(r_cylinder*Theta*sind(beta)+ radius_to_extra.*cosd(beta)-radius_to_extra.*cos(phi))- T(5)*Theta);
Denom_4 = (m_t+(I/r_cylinder))+(m_extra*(1+2*(radius_to_extra^2/r_cylinder^2)*(1+cos(Theta)))+ 2*(radius_to_extra^2/r_cylinder^2));
Quotient_4 = Num_4./Denom_4;
V_c_4= sqrt(Quotient_4);

w = V_c_4/r_cylinder;

end

