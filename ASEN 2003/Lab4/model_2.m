function [w] = model_2(theta)
% Purpose: calculations for model 2 - models a balance cylinder with a support structure as it 
% rolls down a ramp with angle beta and accounts for a negative moment from friction in the system.
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
T = [1 1.2 1.4 1.45 1.5]; % torque friction of bearing
m_t = m_cylinder + m_supports; % [kg]

%% Model 2
% Pre allocate vector
w_2store = [];

% For loop runs through length of 'T', torque friction of bearing
for i = 1:length(T)
    
Num_2 = ((2*m_t*g*r_cylinder*theta*sind(beta))- T(i)*theta);
Denom_2 = (m_t + I/(r_cylinder^2));
Quotient_2 = Num_2./Denom_2;
v_g_2 = sqrt(Quotient_2);

w_2 = v_g_2/r_cylinder;

w_2store = [ w_2store w_2];
end

w = w_2store;

end