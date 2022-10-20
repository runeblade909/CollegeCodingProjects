clc
clear


% Loop track section

% This function takes the xyz position of the start of the track section,
% the incoming velocity vector, the desired radius of the loop, and number
% between 0 and 1 representing what percentage of the loop to be completed
% is desired (so you could do a half loop for instance).  It outputs a
% vector of xyz coordinates for the loop section, the velocity vector
% leaving the loop, and the g forces for the entire loop. 
% The particle is on the inside of the loop

% This function makes a loop in the xz plane
% Theta is defined as zero along the x axis and increasing towards the positive
% z axis
%   z
%   ^
%   |           ^
%   |   (theta) |
%   -------->x

xyz= [20 0 20];

[pos,vel_out,g_forces] = loop_xz(xyz,5,.25);

function [xyz, vel_out, g_forces] = loop_xz(xyz_start, vel_in, r, prop)

% grab all of the starting positions from the vector
xin = xyz_start(1);
yin = xyz_start(2);
zin = xyz_start(3);
n = 100;                    %number of points on this track element


% First, compute the necessary starting angle for the loop (referenced like
% a unit circle) based on the incoming velocity vector
vx = vel_in(1);
vz = vel_in(3);
theta_in = atan(vz/vx);


% create theta vector for the loop and the xyz coordinates.  Make sure the
% loop starts at the beginning x, y, and z coordinates
theta = theta_in + linspace(0,2*pi*prop,n);  %Prop 0-1 how much of a loop
x = r * cos(theta) + xin - r * cos(theta_in); % position in the circle plus the shift
y = yin*ones(1,n);                                              %Since it's in the xz plane the y values don't change
z = r * sin(theta) + zin - r *sin(theta_in);


% compute the ending angle to get the velocity
theta_out = theta(end);                              %angle on the circle leaving the loop
h0 = 125;                                                   %starts at 125 m
vel = sqrt(2*9.81*(h0-z));  
dh =  zin - z;                                        %difference in height from start to end of the loop
mag_v_in = sqrt(vx.^2+vz.^2);

speed =  mag_v_in + sqrt(abs(2*g.*(zin-z)));                                 %speed coming out of the loop
vel = [(speed.*cos(theta))' (0.*cos(theta))' (speed.*sin(theta))'];

vel_out = [vel(end,1) 0 vel(end,3)];    %Velocity vector leaving the loop


% Now compute the g-forces experienced at every point on the loop for the particle on the inside of the loop                %Speed during every point on the loop
g = 9.81;                                   %Force of gravity [m/s^2]
gs = -1 + vel.^2 * (r*g) ;     %G-forces experienced, derivation written up

% Compile the g forces and xyz coordinates into the matrices to be outputted
g_forces = [(gs(:,1).*cost(theta))' (0.*cos(theta))' (gs(:,3).*sin(theta))'];   % G-force matrix [front/back, left/right, up/down]
xyz = [x' y' z'];                           %XYZ matrix

end