function [aircraft_state_trim, control_surfaces_trim] = initialTrim(trim_definition,trim_variables)
%
%
% INPUT:    trim_definition = [V0; h0]
%
%           trim_variables = [alpha0; de0; dt0];
%
% OUTPUT:   aircraft_state_trim = [xe, ye, ze, roll, pitch, yaw, ue, ve, we, p, q, r]
%           control_surfaces_trim = [de da dr dt]
%
% 
% METHOD:   Determines initial state trim conditions and the required
%           control surfaces to stay in trim.

u = trim_definition(1) * cos(trim_variables(1));
w = trim_definition(1) * sin(trim_variables(1));

aircraft_state_trim = [0,0,-trim_definition(2),0,trim_variables(1),0,u,0,w,0,0,0]';

control_surfaces_trim = [trim_variables(2),0,0,trim_variables(3)]';





end

