%% Names: Zak Reichenbach, Madison Ritsch, Jordan Brownlow

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CODE CHALLENGE 1 - 
%
% The purpose of this challenge is to estimate atmospheric pressure in
% Boulder CO using a pressure model and measurements, and compare the two
% through error analysis and statistics. 
%
% To complete the challenge, execute the following steps:
% 1) Load the given dataset
% 2) Extract altitude and pressure data
% 3) Determine standard deviation, variance, mean, and 
%    standard error of the mean of the pressure data
% 4) Using information given about the instrument, find uncertainty associated
%    with altitude measurements
% 5) Use the model to predict pressure measurements at each altitude in the
%    data set, along with propagated uncertainty
% 6) Compare results, discuss, and print answers to the command window. 
% Bonus) Repeat for larger measurement uncertainty in altitude
%
% NOTE: DO NOT change any variable names already present in the code.
% 
% Upload your team's script to Canvas to complete the challenge.
% 
% NAME YOUR FILE AS Challenge1_Sec{section number}_Group {group breakout #}.m 
% ***Section numbers are 1 or 2*** 
% EX File Name: Challenge1_Sec1 _Group15.m 
%
%
% 1) 
% 2) 
% 3) 
% 4) 
% 5)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Housekeeping
clear all   % Clear all variables in workspace
close all   % Close all open figure windows
clc           % Clear the command window

%% 1) Load data from given file 
Data = readtable('PressureInBoulder.csv');

%% 2) Extract just the altitude and station pressure data columns to meaningfully named variables

AltitudeData = table2array(Data(:,3));
h = AltitudeData;
PressureData = table2array(Data(:,2));
% 
%% 3) Determine Statistics and Error
% % the standard deviation, variance, mean, and standard error of the mean (sem) of the pressure data
N = length(PressureData);

MeanPressure = sum(PressureData) / N;

BLAH = sum(PressureData-MeanPressure);
% 
StdevPressure = sqrt((1/(N))*(BLAH^2));
% 
VarPressure = (1/(N))*(BLAH^2);
% 

% 
Sem_Pressure= StdevPressure/sqrt(N); 
% 
%% 4) Uncertainty
% %    The altitude measurements were taken using an instrument that displayed
% %    altitude to the nearest tenth of a meter. 
% 
% %    What is the associated absolute uncertainty with these measurements? 
% 
AltitudeUncertainty = 0.05;   % [m]
% 
%% 5) Pressure Predictions 
% %    Using the altitude measurements and uncertainty, predict pressure with the follwing model:
% %    First, propagate uncertainty BY HAND before calculating uncertainty for each value. 
% %    Then check: is it different for each calculation?
% 
% %               Model
% %    P_est =   P_s * e^(-k*h)               
% %    Assume P_s is 101.7 ± 0.4 kPa and k is 1.2*10^(-4) [1/m] 
% 
 P_s = 101.7;       % ± 0.4 [kPa]
 k = 1.2*10^(-4);   % [1/m]
% 
 P_est = P_s*exp(-k*AltitudeData);
% GENERAL RULE QUADRITURE
 P_sig = sqrt((exp(-k.*h)*(.4)).^2+(P_s.*(-h).*exp(-k.*h)*(5*10^(-6))).^2+(P_s*(-k)*exp(-k.*h)*(AltitudeUncertainty)).^2);

 % 
%% 6) Print Results
% %    Display the predicted pressure from the model with it's associated uncertainty and
% %    the average pressure with the it's standard error of the mean from the data.
% 
 results  = table(P_est,P_sig);
 P_data = [num2str(MeanPressure) '  ± ' num2str(Sem_Pressure) '   kPa'];
 disp(results);
 disp(P_data);
% 
% % Disucss the accuracy of the model and whether or not you think the
% % model agrees with the measurements 
% 
 fprintf('\nModel Discussion: If we look at the printed vector variable "results", \nwe see that the estimated pressure with the added uncertainty comes well within \nour mean pressure and uncertainty. So we collective agree that the estimate \nwe made based on the model agree with the averages.\n ')
% 
%% Bonus
% %   Repeat steps 4-6, but assume the altitude measurements were taken on a
% %   lower precision instrument that only displayed altitude to nearest 10
% %   meters
% %   How does this change the results and comparison ? 
% 
% altitude_uncertainty_new = 5   % [m]

            %% 4) Uncertainty
        % %    The altitude measurements were taken using an instrument that displayed
        % %    altitude to the nearest tenth of a meter. 
        % 
        % %    What is the associated absolute uncertainty with these measurements? 
        NewAltitudeData = round(AltitudeData,-1);
        AltitudeUncertainty = 5;   % [m]
        % 
        %% 5) Pressure Predictions 
        % %    Using the altitude measurements and uncertainty, predict pressure with the follwing model:
        % %    First, propagate uncertainty BY HAND before calculating uncertainty for each value. 
        % %    Then check: is it different for each calculation?
        % 
        % %               Model
        % %    P_est =   P_s * e^(-k*h)               
        % %    Assume P_s is 101.7 ± 0.4 kPa and k is 1.2*10^(-4) [1/m] 
        % 
         P_s = 101.7;       % ± 0.4 [kPa]
         k = 1.2*10^(-4);   % [1/m]
        % 
         P_est = P_s*exp(-k*NewAltitudeData);
        % GENERAL RULE QUADRITURE
         P_sig = sqrt((exp(-k.*h)*(.4)).^2+(P_s.*(-h).*exp(-k.*h)*(5*10^(-6))).^2+(P_s*(-k)*exp(-k.*h)*AltitudeUncertainty).^2);

         % 
        %% 6) Print Results
        % %    Display the predicted pressure from the model with it's associated uncertainty and
        % %    the average pressure with the it's standard error of the mean from the data.
        % 
         results  = table(P_est,P_sig);
         P_data = [num2str(MeanPressure) '  ± ' num2str(Sem_Pressure) '   kPa'];
         disp(results);
         disp(P_data);
        % 
        % % Disucss the accuracy of the model and whether or not you think the
        % % model agrees with the measurements 
        % 
         fprintf('\nModel Discussion BONUS: When we assumed that the altitude measurements were taken on a lower precision instrument that only displayed altitude\n to the nearest 10 meters, our model was still relatively accurate. However, as compared to a more accurate measuring device,\n the results are farther off of the mean for the less accurate device. When thinking intuitively, this would make sense as a more accurate\n device has a smaller standard deviation and a smaller error.')
        % 
