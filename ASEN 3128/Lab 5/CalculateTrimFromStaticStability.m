% Eric W. Frew
% ASEN 3128
% CalculateTrimFromStaticStability.m
% Created: 10/15/20
% STUDENTS COMPLETE THIS FUNCTION

function [alpha_trim, elevator_trim] = CalculateTrimFromStaticStability(trim_definition, aircraft_parameters)
%
% Inputs:	trim_definition         = [V0; h0]
%           aircraft_parameters     = structure with A/C parameters
%
%
% Outputs:	alpha_trim
%           elevator_trim
%
%
% Methodology: Uses static stability analysis to set up matrices and solve
% for trim angle of attack and elevator defelction.


%%% I know, I usually do not like doing this, but it makes the problem
%%% easier to understand in this case
ap = aircraft_parameters; % redefine states and inputs for ease of use
Va_trim = trim_definition(1);
h_trim= trim_definition(2);
rho_trim = stdatmo(h_trim);


%%% Determine lift coefficient needed for trim
CLtrim = ap.W / (0.5 * rho_trim * Va_trim^2 * ap.S); %STUDENTS COMPLETE

%%% Set up matrices and vector relationship

CoefMat = [ap.CLalpha,ap.CLde;ap.Cmalpha,ap.Cmde];
trimvec = [CLtrim;-ap.Cm0];
    
%%% Solve for angle of attack and elevator angle from matrix and vector

result = (CoefMat)\trimvec;
alpha_trim = result(1);
elevator_trim = result(2);