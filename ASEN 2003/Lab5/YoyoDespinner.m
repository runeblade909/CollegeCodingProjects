% ASEN 2003 Yoyo Despinner Exercise]
% Zak Reichenbach
%3/18/2021

%% House Keeping
clc
clear all
close all

%% Data Import

%Data looks like [Time[ms], Speed[rpm]]

datafiles = dir('RopeData\');
num_files = length(datafiles)-2;
data = cell(1,num_files);
names = strings(1,num_files);

for i = 1:num_files
    % Add 3rd col of file data to second col of each cell
    tempdata = readmatrix(strcat('RopeData\', datafiles(i+2).name));
    data{1,i} = tempdata;
    
end

%% Satellite/Mass Properties
Is = 0.0063; % +- 0.0001 kg*m^2
R = 0.076; % M (3 inches)
M = 54; %g
% 1 inch = 0.0254m;

%% Tangential Model

% We want to calculate the length of chord needed for a tangential release




M2 = (2*M)/(1000);
Lt = sqrt(R^2+(Is/M2)); %m;
Lt_Inches = Lt*(100/2.54);

Lt_Actual = Lt_Inches/2;

C = Is/(M2*R^2)+1;


Ltt = R*sqrt(C);
Ltt_Inches = Ltt/0.0254;

% NASA Radial Cord Length

Lr = R*(sqrt(Is/(M2*R^2)+1)-1);
Lr_Inches = Lr/0.0254;


Derived_Lr = (Is+M2*R^2-M2*R*sqrt((Is/M2)+R^2))/(M2*sqrt((Is/M2)+R^2));




%% Plots
plots = 1;
if plots == 1


        hold on
        plot(data{1,1}(:,1),data{1,1}(:,2))
%         plot(data{1,2}(:,1),data{1,2}(:,2))
%         plot(data{1,3}(:,1),data{1,3}(:,2))
        plot(data{1,4}(:,1),data{1,4}(:,2))
        plot(data{1,5}(:,1),data{1,5}(:,2))
        plot(data{1,6}(:,1),data{1,6}(:,2))
%         plot(data{1,7}(:,1),data{1,7}(:,2)) 
        plot(data{1,8}(:,1),data{1,8}(:,2))
        plot(data{1,9}(:,1),data{1,9}(:,2))
        yline(0)
        hold off
        
        title('RPM vs Time')
        xlabel('Time [ms]')
        ylabel('RPM')
        legend('5.5','7','7.5','8','9','No Masses')

end