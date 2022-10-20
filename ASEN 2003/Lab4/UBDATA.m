function [theta_exp, w_exp, time] = UBDATA(filename)
% Purpose: Loads and cleans data for 'UBMAIN.m'
% Authors: Jacob Starkel, Krystal Horton, Cate Billings, and Zak Reichenbach
% Date Completed: 3/12/2021

%Load File
file = load(filename);

% Clean Data by breking file into arrays for time, theta, and w
time = file(:,1);
theta_exp = file(:,2);
w_exp = file(:,3);
end