%% main script

% Zak Reichenbach

clear
clc

% initialize constants
L = 27.25;
w_0 = 2001;

error = zeros(1,30);

% for loop that varies p from 1 to 30
for p = 1:30
    % call discretize_load and pass output to moment_error
    matrix = discretize_load(p, L, w_0);
    
    [r_F, r_M] = wall_reactions(matrix);

    % save error in a vector
    error(p) = abs(moment_error(matrix, L, w_0));
end

% plot p vs error
p = 1:30;
plot(p, error);
xlabel('Number of Point Loads');
ylabel('Error');
title('Number of Point Loads vs. Error');

%% Function that calculates the resultant force of each point load and the 
% distance of the point load from the origin

function[matrix]= discretize_load(p, L, w_0)

right_force = zeros(p,2);
matrix = zeros(p,2);

% Divide L by p so find increment amount
incr = L/p;

x_i = incr;

% Create for loop to store values of x_i and calculate right height of trap
for row = 1:p
    
        matrix(row,2) = x_i;

        % Calculate force of right side
        right_force(row,1) = w_0*(1-(x_i./L));

        x_i = x_i + incr;
        
end

% calculate area of trap to calculate resultant force using geometry
for k = 1:p
    
    if k == 1
        matrix(k,1) = incr * ((w_0 + right_force(k,1))/2);
    else
        matrix(k,1) = incr * ((right_force(k,1) + right_force(k-1,1))/2);
    end
    
end

end


%% Calculates the reactions at the wall given point forces

function[r_F, r_M] = wall_reactions(matrix)
    
    L = 27.25;
    
    F_Ri = sum(matrix(:,1));
    % A_y = (w_0*L)/2 = 2.7264 x 10^4
    
    r_F = sum(F_Ri);
    
    % M_A = (w_0*L^2)/6 = 2.4764 x 10^5
    r_M = (1/3)*L*F_Ri;

end


%% function calculates the bending moment at a randomly-chosen point 
% x = 3L/16 caused by series of point loads and compare to bending 
% moment found with the distributed load from Part 2

function[error] = moment_error(matrix, L, w_0)

x = 3*L/16;

% Finds value closest to x = 3L/16
v = sort(matrix(:));

idx = find(v >= x, 1, 'first');

if ~isempty(idx)
    x_val = v(idx);
end

% finds M_point
M_point = ((L-x_val)/3)*(L-x_val)*(1/2)*(matrix(idx,1));

% M = (w_0/6)*(L-x)^(2)*(1-(x/L)) from part 2

% finds M_dist (simplified with x = 3L/16)
M_dist = (w_0/6)*(13*L/16)^(2)*(13/16);

error = M_point - M_dist;
disp(M_dist);

end