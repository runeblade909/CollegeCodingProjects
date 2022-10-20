%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CODE CHALLENGE 10 - Gaussian Elimination
% 
% This challenge is an exercise in applying Gaussian Elimination in order
% to solve a system of equations. The system of equations you are looking
% to solve is as follows:
%        7x + 3y - 17z = 13
%       -4x + 2z = -2
%        4x + 3y - 9z = -5
%
%
% NOTE: DO NOT change any variable names already present in the code.
% 
% Upload your team's script to Gradescope when complete.
% 
% NAME YOUR FILE AS Challenge10_Sec{section number}_Group{group breakout #}.m 
% ***Section numbers are 1 or 2*** 
% EX File Name: Challenge10_Sec1_Group4.m 
%
% STUDENT TEAMMATES
% 1) Zak Reichenbach
% 2) Anna Casillas
% 3) Jack Iribarren
% 4) Lucas House
% 5) Tristan Workman
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Housekeeping
clear all;
close all;
clc;


%% Organizing Known Values
coeffs1 = [7 3 -17]; % Coefficients of equation 1 in x y z order
coeffs2 = [-4 0 2]; % Coefficients of equation 2 in x y z order
coeffs3 = [4 3 -9]; % Coefficients of equation 3 in x y z order

answers = [13; -2; -5]; % Answers to the three equations

AMat = [coeffs1;coeffs2;coeffs3]; % Creating an A-matrix using the given coefficients
[m,n] = size(AMat); % m for the number of rows in your A-matrix, n for the number of columns



%% Reducing the A-matrix using Recursive Method

%Step 1
%a22 = a22 - (a21/a11)*a12         a32 = a32 - (a31/a11)*a12
%a23 = a23 - (a21/a11)*a13         a33 = a33 - (a31/a11)*a13

%b2 = b2 - (a21/a11)*b1
%b3 = b3 - (a31/a11)*b1

%Step 2 
%a33 = a33-(a32/a22)*a23
%b3 = b3 - (a32/a22)*b2

for j =2:3
   
    for i = 1:3
    AMat(j,i) = 



    end
end


[L,U] = lu(AMat);
% AMat
% L
% U


%% Solving for variables with Back Substitution
...
    
x = nan; % Solving for x
z = nan; % Back-substituting x to solve for z
y = nan; % Back-substituting x and z to solve for y

% Outputting answers
fprintf('The calculated values for our variables are:\n')
fprintf('    x = %.4f\n',x)
fprintf('    y = %.4f\n',y)
fprintf('    z = %.4f\n\n',z)


%% BONUS:
%    Solve the system using matrix methods, with the "\" operator and the
%    "inv()" operator. Compare the time it takes for these two methods with
%    the time that Gaussian elimination took.

% Using \ operator
tic
answersBackslash = nan;
timeBackslash = toc;

% Using matrix inverse
tic
answersInverse = nan;
timeInverse = toc;