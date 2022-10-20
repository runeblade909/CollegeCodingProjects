%% ASEN 3112 Lab 4
% Group: 
% Zach Faith, Zak Reichenbach, Lila Monprode, Joe Miserlian, Gerardo Romero



clc; clear;
load("11_in_rectangle_FIXED.mat");
rect = data;
load("11inch_SquareRod.mat");
square = data; 

%% Experimental Data

% Converting Volts to lbs

conv = (1/2.37)*1000; % conversion factor for load cell

    % Rectangle
    
    V_rec = rect(:,2);
    
    P_rec = conv * V_rec;
    dev = rect(:,3);
    
    plot(dev,P_rec)
    title('Buckling Experiment Rectangular Specimen')
    xlabel('Horizontal Deflection (in)')
    ylabel('Load,P (lbs)')
    
    % Find buckling load:
    idx_cr = rect(:,3) == 1/16;
    Pcr_rec = mean(P_rec(idx_cr));


    % Square
    V_square = square(:,2);
    
    P_square = conv * V_square;
    dev_square = square(:,3);
    
    figure
    plot(dev_square,P_square)
    title('Buckling Experiment Square Rod Specimen')
    xlabel('Horizontal Deflection (in)')
    ylabel('Load,P (lbs)')


    % Find buckling load:
    idx_cr = square(:,3) == 1/16;
    Pcr_square = mean(P_square(idx_cr));
    

    %% Theory


    E = 10e+6; % psi, same for both
    L = 11; % in, same for both
    
    % rectangle
    h_r = 1; % in
    w  = 0.125;% in
    I_r = (h_r*w^3)/12;

    % Buckling Load
    Pcr_rec_theory = pi^2*E*I_r/L^2;
    error_rec = abs(Pcr_rec_theory-Pcr_rec)/Pcr_rec_theory;
    
    disp("Rectangle: ")
    disp("  Experimental P_cr: " + Pcr_rec + " Lbs")
    disp("  Theoretical P_cr: " + Pcr_rec_theory + " Lbs")
    disp("  Error: " + error_rec*100 + "%");
    
    
    
    % square
    side = 0.25; % in
    t = 0.0625; % in
    b =  side  -  2*t;
    I_s = side^4/12 - b^4/12;

    % Buckling Load
    Pcr_square_theory =  pi^2*E*I_s/L^2;
    error_square = abs(Pcr_square - Pcr_square_theory)/Pcr_square_theory;
    
    disp("Square: ")
    disp("  Experimental P_cr: " + Pcr_square + " Lbs")
    disp("  Theoretical P_cr: " + Pcr_square_theory + " Lbs")
    disp("  Error: " + error_square*100 + "%")


