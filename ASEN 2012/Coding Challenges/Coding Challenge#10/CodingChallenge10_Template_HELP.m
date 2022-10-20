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
% NAME YOUR FILE AS Challenge9_Sec{section number}_Group{group breakout #}.m 
% ***Section numbers are 1 or 2*** 
% EX File Name: Challenge9_Sec1_Group15.m 
%
% STUDENT TEAMMATES
% 1) 
% 2) 
% 3) 
% 4)
% 5) 
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Housekeeping
clear all
close all
clc


%% Organizing Known Values
coeffs1 = [nan nan nan]; % Coefficients of equation 1 in x y z order
coeffs2 = [nan nan nan]; % Coefficients of equation 2 in x y z order
coeffs3 = [nan nan nan]; % Coefficients of equation 3 in x y z order

answers = [nan; nan; nan]; % Answers to the three equations

a = []; % Creating an A-matrix using the given coefficients
[m,n] = size(AMat); % m for the number of rows in your A-matrix, n for the number of columns



%% Reducing the A-matrix using Recursive Method
% Initialize the new A matrix, called AA, and b matrix, called bb
aa = a;
bb = answers;

% Using the Recursive Method to reduce our A-matrix, starting at step 1, row 2
r = 1; % step number
for k = 2:m % Total number of iterations that must be taken
    for i = r+1:m % Marker for rows
        for j = r:n % Marker for columns
            scale_factor = a()/a();
            aa(i,j) = a()-a()*scale_factor; % Define the new values of the A matrix
        end
        bb() = answers() - answers()*scale_factor; % Define the new values of the b matrix
    end
    a = aa; % Set the new reference a matrix
    b = bb; % Set the new reference b matrix
    r = r+1; % Increase the step number
end

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