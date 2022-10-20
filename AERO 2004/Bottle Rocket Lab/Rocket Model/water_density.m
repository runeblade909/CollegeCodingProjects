function [density] = water_density(temp, unit)
% WATER_DENSITY
%   DENSITY = water_density(TEMP, UNIT)
%
%   TEMP is type double representing the temperature of water
%
%   UNIT is type char representing the units of temperature TEMP. 
%   Must be 'K', 'C', 'R', or 'F'
%
%   DENSITY is in kg/m^3
%
%   Valid for 0.1 to 373.946 deg C

% arguments
%    temp double
%    unit char
% end
    
% Function to convert temp of water to density.
% Written by Bennett Grow 4/15/2021 for CU ASEN 2004

% Data pulled from Engineering Toolbox 4/15/2021
% https://www.engineeringtoolbox.com/water-density-specific-weight-d_595.html

% ----------------------------------------------------------------------- %

% Check if valid unit and convert to deg C
if (unit == 'K' || unit == 'k')
    temp = temp - 273.15;
elseif (unit == 'C' || unit == 'c')
    % Correct unit
elseif (unit == 'F' || unit == 'f')
    temp = (temp - 32) * (5/9);
elseif (unit == 'R' || unit == 'r')
    temp = (temp - 491.67) * (5/9);
else
    error('Unit must be char: K, C, F, or R.');
end


% Temp data in deg C
temp_table = [0.1; 1; 4; 10; 15; 20; 25; 30; 35; 40; 45; 50; 55; 60; 65; 70; 75; 80; 85; 90; 95; 100; 110; 120; 140; 160; 180; 200; 220; 240; 260; 280; 300; 320; 340; 360; 373.946];
% Density data in g/cm^3
density_table = [0.9998495; 0.9999017; 0.9999749; 0.9997000; 0.9991026; 0.9982067; 0.9970470; 0.9956488; 0.9940326; 0.9922152; 0.99021; 0.98804; 0.98569; 0.98320; 0.98055; 0.97776; 0.97484; 0.97179; 0.96861; 0.96531; 0.96189; 0.95835; 0.95095; 0.94311; 0.92613; 0.90745; 0.88700; 0.86466; 0.84022; 0.81337; 0.78363; 0.75028; 0.71214; 0.66709; 0.61067; 0.52759; 0.3220];
% Convert density data to kg/m^3
density_table = density_table * 1000;

if (temp < temp_table(1))
    warning('Model only valid for 0.1 to 373.946 degC. Returning closest density.');
    density = density_table(1);
elseif (temp > density_table(end))
    warning('Model only valid for 0.1 to 373.946 degC. Returning closest density.');
    density = temp_table(end);
else
    % Linerally interpolate between Enginering Toolbox data
    density = interp1(temp_table, density_table, temp);
end

end