
%% Appendix2DataProcessing.m
%{

Created by Bennett Grow
ASEN 2004 Lab Milestone 2, Spring 2021, CU Boulder

Reads in data from appendix 2 graphs, polyfits the data, then saves a
single matrix app2data

app2data = [AOA CL CD Cm4]

%}

%% Read in raw graph data
CL_AOA_raw = readmatrix("Appendix2Data.xlsx", "sheet", 1);
CD_AOA_raw = readmatrix("Appendix2Data.xlsx", "sheet", 2);
Cm4_AOA_raw = readmatrix("Appendix2Data.xlsx", "sheet", 3);

%% Polyfit all data
AOA = -20:30;
fit = 6;

% CL vs AOA
pCL = polyfit(CL_AOA_raw(:,1), CL_AOA_raw(:,2), fit);
CL = polyval(pCL, AOA);

% CD vs AOA
pCD = polyfit(CD_AOA_raw(:,1), CD_AOA_raw(:,2), fit);
CD = polyval(pCD, AOA);

% Cm4 vs AOA
pCm4 = polyfit(Cm4_AOA_raw(:,1), Cm4_AOA_raw(:,2), fit);
Cm4 = polyval(pCm4, AOA);

% Rotate data
AOA = AOA';
CL = CL';
CD = CD';
Cm4 = Cm4';

%% Plot
plotgraphs = false;

if plotgraphs == true
    
    % CL
    figure
    scatter(CL_AOA_raw(:,1), CL_AOA_raw(:,2));
    hold on
    grid on
    plot(AOA, CL);
    title("CL vs AOA")
    hold off
    
    % CD
    figure
    scatter(CD_AOA_raw(:,1), CD_AOA_raw(:,2));
    hold on
    grid on
    plot(AOA, CD);
    title("CD vs AOA")
    hold off
    
    % Cm4
    figure
    scatter(Cm4_AOA_raw(:,1), Cm4_AOA_raw(:,2));
    hold on
    grid on
    plot(AOA, Cm4);
    title("Cm4 vs AOA")
    hold off


end


%% Save all data
% Most accurate fit in range -15:20 deg AOA (rows 6 - 41)
range = 6:41;

app2data = [AOA(range) CL(range) CD(range) Cm4(range)];

save('app2data', 'app2data');
clear all;


