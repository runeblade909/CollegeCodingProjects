function [w] = model_1(theta)
% Purpose: calculations for model 1 - models a balance cylinder with a support structure as it
% rolls down a ramp with angle beta
% Authors: Jacob Starkel, Krystal Horton, Cate Billings, and Zak Reichenbach
% Date Completed: 3/12/2021

%% Givens 
m_cylinder  = 11.7; % [kg]
m_supports  = 0.7;  % [kg]
r_cylinder  = .235;  % [m]
k = .203; % [m]
I = m_cylinder*k^2; % [kgm^2]
beta = 5.5; % slope of ramp [degrees]
g = 9.81; % [m/s^2]
m_t = m_cylinder + m_supports; % [kg]


%% Model 1
Num = (2*m_t*g*r_cylinder*theta*sind(beta));
Denom = (m_t + I/(r_cylinder^2));
Quotient = Num/Denom;
v_g = sqrt(Quotient);

w = v_g/r_cylinder;

end