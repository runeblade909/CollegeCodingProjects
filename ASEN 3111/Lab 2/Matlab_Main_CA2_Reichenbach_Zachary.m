%% ASEN 3111 - Computational Assignment 2 - Main
% Model the flow characteristics over an airfoil with vortex panel method 
% along a given chord length by finding and calculating the stream 
% function, velocity potential, and pressure contour

% Author: Zak Reichenbach
% Collaborators: Ketan Kamat
% Date: 10/11/2021

%% Housekeeping

clc;
clear;
close all;
tic

%% Constant Declaration

c = 1.5;            % chord length [m]
alpha = deg2rad(6); % angle of attack [degrees]
V_inf = 30;         % free stream velocity [m/s]
P_inf = 101.3*10^3; % free stream pressure [Pa]
rho_inf = 1.225;    % free stream density [kg/m^3]
N = 100;           % number of vortices

%% Call function to find S.F., V.P., and P.C. for the desired flow
% last input determines what runs in the function:
%   1 - regular function that calculates S.F., V.P., and P.C. for 50000
%       vortices and then finds the number of vortices to get less 
%       than a 0.5% error (Bullet point 2)
%   2 - regular function that is used to conducts study to determine 
%       changes in angle of attack


%% Problem 1



%Plot the base case with N = 25
[x0,y0,Psi0,Phi0,P0,levels0,bounds0,~,~] = Plot_Airfoil_Flow(c,alpha,V_inf,P_inf,rho_inf,10,3);

% % %Look for error vs. model of assumed true value (N = 50,000)
[x,y,Psi,Phi,P,levels,bounds,N2,err] = Plot_Airfoil_Flow(c,alpha,V_inf,P_inf,rho_inf,N,1);

% % %Plot for high angle of attack with more accurate case (2600)
[x1,y1,Psi1,Phi1,P1,levels1,bounds1,~,~] = Plot_Airfoil_Flow(c,deg2rad(30),V_inf,P_inf,rho_inf,N2,3);

    %%%%%%%%%%%%%%% To save on run time
    filename = 'FlatPlate_CA2_Zak_Reichenbach';
    save(filename,'x','y','Psi','Phi','P','levels','bounds','N2','err','x1','y1','Psi1','Phi1','P1','levels1','bounds1');
    load('FlatPlate_CA2_Zak_Reichenbach');

    
%     %% Plots, since other plots in the function are suppressed for run time
%     
%     %% Plot 2
%     figure
%     contourf(x,y,Psi,75) % stream lines
%     hold on 
%     plot([0 c],[0 0],'k','linewidth',2) % airfoil
%     axis equal
%     axis(bounds)
%     ylabel('y')
%     xlabel('x')
%     title(['Stream Lines for N = ' num2str(N2) ' Vorticies at \alpha = 6^{o}']);
% 
%     % Plot equipotential lines at levels
%     figure
%     contourf(x,y,Phi,100) % equipotential lines
%     hold on 
%     plot([0 1.5],[0 0],'k','linewidth',2) % airfoil
%     axis equal
%     axis(bounds)
%     ylabel('y')
%     xlabel('x')
%     title(['Equipotential Lines for N = ' num2str(N2) ' Vorticies at \alpha = 6^{o}']);
% 
%     % Plot pressure contours
%     figure
%     contourf(x,y, P,100) % pressure contour
%     hold on
%     plot([0 c],[0 0],'k','linewidth',2) % airfoil
%     axis equal
%     axis(bounds)
%     k = colorbar;
%     k.Label.String = 'Pressure [Pa]';
%     ylabel('y')
%     xlabel('x')
%     title(['Pressure Contours for N = ' num2str(N2) ' Vorticies at \alpha = 6^{o}']);
%         
%     
%      %% Plot 3
%     figure
%     contourf(x1,y1,Psi1,75) % stream lines
%     hold on 
%     plot([0 c],[0 0],'k','linewidth',2) % airfoil
%     axis equal
%     axis(bounds)
%     ylabel('y')
%     xlabel('x')
%     title(['Stream Lines for N = ' num2str(N2) ' Vorticies at \alpha = 30^{o}']);
% 
%     % Plot equipotential lines at levels
%     figure
%     contourf(x1,y1,Phi1,100) % equipotential lines
%     hold on 
%     plot([0 1.5],[0 0],'k','linewidth',2) % airfoil
%     axis equal
%     axis(bounds)
%     ylabel('y')
%     xlabel('x')
%     title(['Equipotential Lines for N = ' num2str(N2) ' Vorticies at \alpha = 30^{o}']);
% 
%     % Plot pressure contours
%     figure
%     contourf(x1,y1, P1,100) % pressure contour
%     hold on
%     plot([0 c],[0 0],'k','linewidth',2) % airfoil
%     axis equal
%     axis(bounds)
%     k = colorbar;
%     k.Label.String = 'Pressure [Pa]';
%     ylabel('y')
%     xlabel('x')
%     title(['Pressure Contours for N = ' num2str(N2) ' Vorticies at \alpha = 30^{o}']);
%         
%     

    % Plot error
    N_spacing = 100:100:N2;
    figure
    plot(N_spacing,err(1:length(N_spacing))*100,'linewidth',2);
    ylabel('Percent Error')
    xlabel('Number of Vortices')
    title(['Pressure Contours Error for N = ' num2str(N2) ' Vorticies']);

fprintf('As seen in figures 1-3 and the higher resolution figures 4-6, more panels produce \nplots that are not only closer to reality, but make more sense with lines \nthat we would expect around our airfoils with more accurate values and straighter lines.\n\n')
fprintf('Figure 10 displays how the component with most error out of the stream function, \nvelocity potential, and pressure contour, and how it changes as the number of\nvortices increases. It was found that %i is the number of\nvortices required to have a error less than 0.5%%. This\naccuracy was based off of a 50,000 vortice model being taken as the correct\ntheoretical value and error was determined from there. In the end,\nas the number of votices increases the better the model becomes at\nmodelling the system. Theoretically the best model would have an\ninfinite amount of vorticies but that is not feasible.\n\n',N2);
fprintf('As seen in our pressure contor figures 3 & 6, it is apparent that we have the highest \npressure on the side of the wing directly facing the uniform flow with the \nlowest pressure being on the other side of the airfoil. With increased angle of \nattack as seen in figures 7-9, we can see some seperation in our stream and \nequipotential lines as well as a general large increase in pressure.\n\n')
%% Study on how angle of attack changes S.F., V.P., and P.C.
% last input determines what runs in the function:
%   1 - regular function that calculates S.F., V.P., and P.C. for 50000
%           vortices and then finds the number of vortices to get less 
%           than a 0.5% error
%   2 - regular function that is used to conducts study to determine 
%           changes in 1) chord length, 2) angle of attack, and 3) free
%   3 - regular function that calculates S.F., V.P., and P.C. for the
%           N input value of vortices
%
%
% angle of attack of -6 deg
alpha_a_Neg6 = deg2rad(-6);
alpha_a_0 = 0;
alpha_a_6 = deg2rad(6);
alpha_a_9 = deg2rad(9);
[x_a_Neg6,y_a_Neg6,Psi_a_Neg6,Phi_a_Neg6,~,~,bounds_a_Neg6,~] = Plot_Airfoil_Flow(c,alpha_a_Neg6,V_inf,P_inf,rho_inf,N,2);

% angle of attack of 6 deg
[x_a6,y_a6,Psi_a6,Phi_a6,~,~,bounds_a6,~] = Plot_Airfoil_Flow(c,alpha_a_6,V_inf,P_inf,rho_inf,N,2);

% angle of attack of 9 deg
alpha_a9 = deg2rad(9);
[x_a9,y_a9,Psi_a9,Phi_a9,~,~,bounds_a9,~] = Plot_Airfoil_Flow(c,alpha_a9,V_inf,P_inf,rho_inf,N,2);

% plot S.F.
figure
subplot(3,1,1)
contourf(x_a_Neg6,y_a_Neg6,Psi_a_Neg6,40) % stream function
hold on 
plot([0 c],[0 0],'k','linewidth',2) % airfoil
axis(bounds_a_Neg6)
title(['Stream Lines for N = ' num2str(N) ' Vorticies at \alpha = ' num2str(-6) '^{o}']);
subplot(3,1,2)
contourf(x_a6,y_a6,Psi_a6,40) % stream function
hold on
plot([0 c],[0 0],'k','linewidth',2) % airfoil
axis(bounds_a6)
title(['Stream Lines for N = ' num2str(N) ' Vorticies at \alpha = ' num2str(6) '^{o}']);
ylabel('y')
subplot(3,1,3)
contourf(x_a9,y_a9,Psi_a9,40) % stream function
hold on
plot([0 c],[0 0],'k','linewidth',2) % airfoil
axis(bounds_a9)
title(['Stream Lines for N = ' num2str(N) ' Vorticies at \alpha = ' num2str(9) '^{o}']);
xlabel('x')
hold off

% plot V.P.
figure
subplot(3,1,1)
contourf(x_a_Neg6,y_a_Neg6,Phi_a_Neg6,100) % velocity potential
hold on 
plot([0 c],[0 0],'k','linewidth',2) % airfoil
axis(bounds_a_Neg6)
title(['Equipotential Lines for N = ' num2str(N) ' Vorticies at \alpha = ' num2str(-6) '^{o}']);
subplot(3,1,2)
contourf(x_a6,y_a6,Phi_a6,100) % velocity potential
hold on
plot([0 c],[0 0],'k','linewidth',2) % airfoil
axis(bounds_a6)
title(['Equipotential Lines for N = ' num2str(N) ' Vorticies at \alpha = ' num2str(6) '^{o}']);
ylabel('y')
subplot(3,1,3)
contourf(x_a9,y_a9,Phi_a9,100) % velocity potential
hold on
plot([0 c],[0 0],'k','linewidth',2) % airfoil
axis(bounds_a9)
title(['Equipotential Lines for N = ' num2str(N) ' Vorticies at \alpha = ' num2str(9) '^{o}']);
xlabel('x')
hold off

fprintf('figures 11 & 12 display how the stream function and the velocity\npotential change as the angle of attack changes. When looking\nat the figure, it is seen that the angle of attack does have an effect\non the stream function or the velocity potential.\n\n');


%% Problem 2
% define the NACA airfoil desired to model
NACA_num = '0012';

% calculate important airfoil values that can be found form the NACA 
% airfoil number
m = str2double(NACA_num(1))/100; % maximum camber
p = str2double(NACA_num(2))/10; % location of maximum camber
t = str2double(NACA_num(3:4))/100; % thickness

% define other characteristics of the airfoil and the flow around it
c = 1; % chord length [m]
V_inf = 30; % free-stream velocity [m/s]
alpha = 5; % angle of attack [deg]

fprintf('See line 157-167 for justification reasoning for error analysis on \nNACA 0012 airfoil for optimum number of vortices. \n\n')

%{ 
 Find the number of panels required to get the coefficient of pressure
 within 1% difference of the previous number of panels. This will be
 done by running a while loop to find the sectional coefficient of lift and
 the coefficient of pressure at a specified number of panels and then
 finding sectional coefficient of lift and the coefficient of pressure
 again for a greater number of panels. Once this is done, the percent
 difference between the two will be found and compared to the desired
 error. If the error found is greater than the error desired, the while loop
 will run again with an increasing number of panels. If the error found is
 less than the desired error, that number of panels will be taken as
 accurate and be used for the rest of the analysis.
%}

% initialization
err = 1;    % initialize error for while loop
N_0 = 10;   % intialize number of employed panels to model the airfoil
N = N_0;
i = 1;      % initialize the vector index


fprintf('Preforming panel analysis for vortex panel method...\n\n');



% % Assumes that 2000 Panels is enough to get a true reading of CL
% [x_Vort,y_Vort] = NACA_Airfoils(m,p,t,c,2000);
% [Cl_Vort, Cp_Vort] = Vortex_Panel(x_Vort,y_Vort,V_inf,alpha,5);

%     %%%%%%%%%%%%%%% To save on run time
%     filename = 'CP_CL_Vortex_CA2_Zak_Reichenbach';
%     save(filename,'Cl_Vort','Cp_Vort','x_Vort','y_Vort');
load('CP_CL_Vortex_CA2_Zak_Reichenbach');
fprintf('2000 Panel Calculation Complete for "truth" data. \n\n');

GoalErr = 0.01;
while err > GoalErr   % 1% error
    
    % increase the number of panels
    N_vec(i) = N;
    
    fprintf('Error iteration vortices: \n%f\n', (N-1)*2);

    
    % find the x and y values of the surface of the airfoil for the second
    % number of panels
    [x_next,y_next] = NACA_Airfoils(m,p,t,c,N); % (x,y) coordinates of the airfoil surface
    
    N = N + N_0/2;
    
    % find sectional coefficient of lift and coefficient of pressure for
    % the second number of panels
    [Cl_next, Cp_next] = Vortex_Panel(x_next,y_next,V_inf,alpha,5); % c_l and C_p of the airfoil
    

    
    %{
    Upon investigating the values of Cp and Cl, the value that improves the
    most with additional vortices is the Cp. The Cl value unreliably goes
    all over the place and is therefore an unreliable variable to determine
    the accuracy of additional vortices. So, to judge the over all
    accuracy of Cp with the change in vortices, the mean was taken between
    a assumed highly accurate case using 2000 vortices, and will be
    compared to a Cp value of increasing accuracy as panels are added from
    10 initial panels. The error difference will then be analyzed over the
    iterations.
    %}
      
    % calculate error
     
      MainDiff = mean(Cp_Vort);
      DiffDiff = mean(Cp_next);
    
%     StartCP(i) = C_p_next(1);
    
    err = abs(abs(MainDiff-DiffDiff)/(MainDiff));
    err_vec(i) = err;
    
    if err < GoalErr
        Cp = Cp_next;
        Cl = Cl_next;
        x = x_next;
        y = y_next;
        
    end
    
    i = i + 1; % increment the vector index
    
end




% print the number of panels needed to the comand window
fprintf('For a %0.2f%% relative error %i panels are needed to model the airfoil\n\n',err*100,(N-1)*2);

% plot airfoil for ideal number of panels calculated
figure
plot(x,y,'linewidth',2);
axis equal
xlabel('x');
ylabel('y');
title(['NACA ' NACA_num ' with ' num2str((N-1)*2) ' panels']);
grid on

%plot error vs number of panels
figure
plot(2*(N_vec-1),err_vec*100,'linewidth',2);
xlabel('Number of Panels');
ylabel('Relative Error [%]');
title(['Relative Error vs Number of Panels (NACA ' NACA_num ') [#2]']);
ylim([0 100])
grid on

%% Problem 2 Angle of Attack Variation

fprintf('Performing angle of attack variation analysis on the coefficinet of pressure and coefficient of lift.\n\n');

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of -6 degrees
[Cl_neg6, Cp_neg6] = Vortex_Panel(x,y,V_inf,-6,1); % c_l and C_p of the airfoil

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of 6 degrees
[Cl_6, Cp_6] = Vortex_Panel(x,y,V_inf,6,1); % c_l and C_p of the airfoil

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of 9 degrees
[Cl_9, Cp_9] = Vortex_Panel(x,y,V_inf,9,1); % c_l and C_p of the airfoil

% plot coefficient of pressure vs angle of attack
figure
subplot(2,2,1)
plot(x(1:end-1)./c,Cp_neg6,'linewidth',2); % -6 deg AoA
set(gca, 'YDir','reverse')
xlabel('x/c');
ylabel('C_{P}');
title('AoA = -6');
xlim([0 1]);
hold on
grid on
subplot(2,2,2)
plot(x(1:end-1)./c,Cp,'linewidth',2); % 0 deg AoA
set(gca, 'YDir','reverse')
xlabel('x/c');
ylabel('C_{P}');
title('AoA = 0');
xlim([0 1]);
hold on
grid on
subplot(2,2,3)
plot(x(1:end-1)./c,Cp_6,'linewidth',2); % 6 deg AoA
set(gca, 'YDir','reverse')
xlabel('x/c');
ylabel('C_{P}');
title('AoA = 6');
xlim([0 1]);
hold on
grid on
subplot(2,2,4)
plot(x(1:end-1)./c,Cp_9,'linewidth',2); % 9 deg AoA
set(gca, 'YDir','reverse')
xlabel('x/c');
ylabel('C_{P}');
title('AoA = 9');
xlim([0 1]);
hold on
grid on
sgtitle(['C_{p} vs Percent Chord (NACA ' NACA_num ')  [#2]']); % title the whole subplot
hold off
% 
figure
plot(x(1:end-1)./c,Cp_neg6,'linewidth',2,'DisplayName','AoA = -6^{o}'); % -6 deg AoA
set(gca, 'YDir','reverse')
xlabel('x/c');
ylabel('C_{P}');
title(['C_{p} vs Percent Chord (NACA ' NACA_num ')  [#2]']);
xlim([0 1]);
hold on
plot(x(1:end-1)./c,Cp,'linewidth',2,'DisplayName','AoA = 0^{o}');
hold on
plot(x(1:end-1)./c,Cp_6,'linewidth',2,'DisplayName','AoA = 6^{o}');
hold on
plot(x(1:end-1)./c,Cp_9,'linewidth',2,'DisplayName','AoA = 9^{o}');
legend
grid on
% 
% plot sectional coefficient of lift vs angle of attack
figure
plot([-6, 0, 6, 9],[Cl_neg6, Cl, Cl_6, Cl_9],'linewidth',2);
xlabel('AoA [deg]');
ylabel('C_{l}');
title(['Sectional Coefficient of lift vs Angle of Attack (NACA ' NACA_num ') [#2]']);
grid on

fprintf('Using the error analysis outline in theory with lines 157-167 in the code, figure 14 shows how error \n coverges over the various iterations. \n\n')

fprintf('As seen in coefficient of pressure figure 18, increasing angle of attack increases the coefficient \nof pressure on the side of the airfoil facing the uniform flow as it did for our flat plate scenario. \n\n')

%% Problem 3

fprintf('Preforming lift slope determination calculations.\n\n');


[Cl, Cp] = Vortex_Panel(x,y,V_inf,0,1);

% Caclulate lift slope
m_L = (Cl_9-Cl_neg6)/(deg2rad(9)-deg2rad(-6)); % lift slope

% Calculate zero-lift angle of attack
AoA_L0 = -(Cl/m_L);

% print lift slope and zero-lift angle of attack values
fprintf('For the NACA %s airfoil\n\tLift slope: %0.3f [1/rad]\n\tZero-lift AoA: %0.3f [deg]\n\n',NACA_num,m_L,rad2deg(AoA_L0));

% plot NACA 0012 sectional coefficient of lift vs angle of attack
figure(34)
plot([-6, 0, 6, 9],[Cl_neg6, Cl, Cl_6, Cl_9],'linewidth',2,'DisplayName',['NACA ' NACA_num]);
hold on
yline(0);
xline(0);

% redefine NACA airfoil
NACA_num = '2412';
% calculate important airfoil values that can be found form the NACA 
% airfoil number
m = str2double(NACA_num(1))/100; % maximum camber
p = str2double(NACA_num(2))/10; % location of maximum camber
t = str2double(NACA_num(3:4))/100; % thickness

% find the x and y values of the surface of the airfoil
[x,y] = NACA_Airfoils(m,p,t,c,N); % (x,y) coordinates of the airfoil surface

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of -6 degrees
[Cl_neg6, ~] = Vortex_Panel(x,y,V_inf,-6,1); % c_l and C_p of the airfoil

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of 0 degrees
[Cl, ~] = Vortex_Panel(x,y,V_inf,0,1); % c_l and C_p of the airfoil

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of 6 degrees
[Cl_6, ~] = Vortex_Panel(x,y,V_inf,6,1); % c_l and C_p of the airfoil

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of 9 degrees
[Cl_9, ~] = Vortex_Panel(x,y,V_inf,9,1); % c_l and C_p of the airfoil

% Caclulate lift slope
m_L = (Cl_9-Cl_neg6)/(deg2rad(9)-deg2rad(-6)); % lift slope

% Calculate zero-lift angle of attack
AoA_L0 = -(Cl/m_L);

% print lift slope and zero-lift angle of attack values
fprintf('For the NACA %s airfoil\n\tLift slope: %0.3f [1/rad]\n\tZero-lift AoA: %0.3f [deg]\n\n',NACA_num,m_L,rad2deg(AoA_L0));

% plot NACA 2412 sectional coefficient of lift vs angle of attack
figure(34)
plot([-6, 0, 6, 9],[Cl_neg6, Cl, Cl_6, Cl_9],'linewidth',2,'DisplayName',['NACA ' NACA_num]);
hold on

% redefine NACA airfoil
NACA_num = '4412';
% calculate important airfoil values that can be found form the NACA 
% airfoil number
m = str2double(NACA_num(1))/100; % maximum camber
p = str2double(NACA_num(2))/10; % location of maximum camber
t = str2double(NACA_num(3:4))/100; % thickness

% find the x and y values of the surface of the airfoil
[x,y] = NACA_Airfoils(m,p,t,c,N); % (x,y) coordinates of the airfoil surface

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of -6 degrees
[Cl_neg6, ~] = Vortex_Panel(x,y,V_inf,-6,1); % c_l and C_p of the airfoil

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of 0 degrees
[Cl, ~] = Vortex_Panel(x,y,V_inf,0,1); % c_l and C_p of the airfoil

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of 6 degrees
[Cl_6, ~] = Vortex_Panel(x,y,V_inf,6,1); % c_l and C_p of the airfoil

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of 9 degrees
[Cl_9, ~] = Vortex_Panel(x,y,V_inf,9,1); % c_l and C_p of the airfoil

% Caclulate lift slope
m_L = (Cl_9-Cl_neg6)/(deg2rad(9)-deg2rad(-6)); % lift slope

% Calculate zero-lift angle of attack
AoA_L0 = -(Cl/m_L);

% print lift slope and zero-lift angle of attack values
fprintf('For the NACA %s airfoil\n\tLift slope: %0.3f [1/rad]\n\tZero-lift AoA: %0.3f [deg]\n\n',NACA_num,m_L,rad2deg(AoA_L0));

% plot NACA 4412 sectional coefficient of lift vs angle of attack
figure(34)
plot([-6, 0, 6, 9],[Cl_neg6, Cl, Cl_6, Cl_9],'linewidth',2,'DisplayName',['NACA ' NACA_num]);
hold on

% redefine NACA airfoil
NACA_num = '2424';
% calculate important airfoil values that can be found form the NACA 
% airfoil number
m = str2double(NACA_num(1))/100; % maximum camber
p = str2double(NACA_num(2))/10; % location of maximum camber
t = str2double(NACA_num(3:4))/100; % thickness

% find the x and y values of the surface of the airfoil
[x,y] = NACA_Airfoils(m,p,t,c,N); % (x,y) coordinates of the airfoil surface

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of -6 degrees
[Cl_neg6, ~] = Vortex_Panel(x,y,V_inf,-6,1); % c_l and C_p of the airfoil

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of 0 degrees
[Cl, ~] = Vortex_Panel(x,y,V_inf,0,1); % c_l and C_p of the airfoil

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of 6 degrees
[Cl_6, ~] = Vortex_Panel(x,y,V_inf,6,1); % c_l and C_p of the airfoil

% find sectional coefficient of lift and coefficient of pressure for an
% angle of attack of 9 degrees
[Cl_9, ~] = Vortex_Panel(x,y,V_inf,9,1); % c_l and C_p of the airfoil

% Caclulate lift slope
m_L = (Cl_9-Cl_neg6)/(deg2rad(9)-deg2rad(-6)); % lift slope

% Calculate zero-lift angle of attack
AoA_L0 = -(Cl/m_L);

% print lift slope and zero-lift angle of attack values
fprintf('For the NACA %s airfoil\n\tLift slope: %0.3f [1/rad]\n\tZero-lift AoA: %0.3f [deg]\n\n',NACA_num,m_L,rad2deg(AoA_L0));

% % plot NACA 2424 sectional coefficient of lift vs angle of attack
figure(34)
plot([-6, 0, 6, 9],[Cl_neg6, Cl, Cl_6, Cl_9],'linewidth',2,'DisplayName',['NACA ' NACA_num]);
xlabel('AoA [deg]');
ylabel('C_{l}');
title('NACA Airfoils Lift Slopes [#3]');
hold on
grid on
legend('location','northwest');
hold off

fprintf('As seen in figure 34 (NOT AT THE BOTTOM OF THE LIST OF FIGURES), \nthe larger value of camber, the more negative \nthe zero lift angle of attack is. The thickness of a airfoil on the \nother hand gives a steaper lift slope. \n\n')

fprintf('For clarification sake, code ran on the authors machine in 56.455284 seconds. \n\n')

toc