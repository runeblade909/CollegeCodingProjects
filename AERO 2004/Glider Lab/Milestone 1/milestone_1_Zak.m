clc
clear

%% Data Collection:
TempestData = ...
    readmatrix('Tempest UAS & B747 Airfoil and CFD Data for ASEN 2004 Aero Lab (UPDATED).xlsx','Sheet',2);
TempestRaw = ...
    readmatrix('Tempest UAS & B747 Airfoil and CFD Data for ASEN 2004 Aero Lab (UPDATED).xlsx','Sheet',3);
BoeingData = ...
    readmatrix('Tempest UAS & B747 Airfoil and CFD Data for ASEN 2004 Aero Lab (UPDATED).xlsx','Sheet',5);
BoeingRaw = ...
    readmatrix('Tempest UAS & B747 Airfoil and CFD Data for ASEN 2004 Aero Lab (UPDATED).xlsx','Sheet',6);
%% Whole Aircraft Drag Polar Tempest Raw Data
TAoA_Raw = TempestRaw(:,1);
TCL_Raw = TempestRaw(:,2);
TCD_Raw = TempestRaw(:,3);


%% Wholare Aircraft Drag Polar Boeing Raw Data
BCL_Raw = BoeingRaw(:,13);
BCD_Raw = BoeingRaw(:,14);


%% Constants
%Span Efficiency
e = 0.9;

%% Tempest 2D and 3D Lift Curves

%2D Lift Curve Data
TAoA = TempestData(1:18,1);
T_Cl = TempestData(1:18,2);
T_Cd = TempestData(1:18,3);
T_AR = 16.5;

%Creating a 3D Lift Curve

%Creating a Line Fit
[p,~] = polyfit(TAoA(1:14),T_Cl(1:14),1);
Ta0 = p(1);
Ta = Ta0/(1+((57.3*Ta0)/(pi*e*T_AR)));
T_CL = Ta * (TAoA(1:14) - TAoA(4));

%Wing Drag Polar equation(3)
T_CD = T_Cd(1:14) + (T_CL).^2/(pi*e*T_AR);

%Plotting Graphs
%Graph: Angle of Attack against 2D and 3D Lift Curve for Tempest
figure(1)
hold on
grid on
plot(TAoA,T_Cl,'LineWidth',2)
plot(TAoA(1:14),T_CL,'LineWidth',2)
xlabel('Angle of Attack (Deg)')
ylabel('Coefficient of Lift')
title('Tempest 2D and 3D Lift Curve')
legend('2D Lift Curve','3D Lift Curve','Location','SouthEast')
hold off

%% Boeing 2D and 3D Lift Curves
BAoA = TempestData(1:19,1);
B_Cl = BoeingData(1:19,2);
B_Cd = BoeingData(1:19,3);
B_AR = 7;

%Creating a 3D Lift Curve

%Creating a Line Fit
[p1,~] = polyfit(BAoA(1:15),B_Cl(1:15),1);
Ba0 = p1(1);
Ba = Ba0/(1+((57.3*Ba0)/(pi*e*B_AR)));
B_CL = Ba * (BAoA(1:15) - BAoA(4)); 

%Wing Drag Polar equation(3)
B_CD = B_Cd(1:15) + (B_CL).^2/(pi*e*B_AR);

%Plotting Graphs
%Graph: Angle of Attack against 2D and 3D Lift Curve for Boeing
figure(2)
hold on
grid on
plot(BAoA,B_Cl,'LineWidth',2)
plot(BAoA(1:15),B_CL,'LineWidth',2)
xlabel('Angle of Attack (Deg)')
ylabel('Coefficient of Lift')
title('Boeing 2D and 3D Lift Curve')
legend('2D Lift Curve','3D Lift Curve','Location','SouthEast')
hold off

%% Tempest Whole Aircraft Drag Polar
%Skin Friction Coefficient (Glider/Civil Transport)
TCfe = 0.0030;
%Reference Area
TSref = .668; 
%Wetted Area
TSwet = 2.4;  
%Calculating CDmin for Tempest (10)
TCDmin = TCfe*(TSwet/TSref); 
%Calculating Oswald's efficiency for Tempest (12)
Te0 = 1.78*(1-0.045*T_AR^0.68)-0.64; 
%Calculating k for Tempest  %(5)
Tk1 = 1/(pi*Te0*T_AR);
%Going through the 2D coefficient of drag data to find the closest min
%value
[~,i] = min(abs(TCDmin-T_Cd));
%Finding the corresponding CLmin
TCLmin = T_CL(i);

Tk2 = -2 * Tk1 * TCLmin; %Check out later

%Parasite Drag Calculation for Tempest
TCD0 = TCDmin + (Tk1 * TCLmin^2); 
%Whole Aircraft Drag Polar for Tempest
TCD = (TCD0 + (Tk1 * T_CL.^2));

%Plotting Graphs
%Graphs: Whole Aircraft Drag Polar (True and Model) along with Wing Drag
%Polar
figure
hold on
grid on
plot(T_CL,T_CD,'LineWidth',2)
plot(TCL_Raw,TCD_Raw,'LineWidth',2)
plot(T_CL,TCD,'LineWidth',2)
hold off
legend('Approximated Wing Drag Polar',' True Whole Aircraft Drag Polar','Approximated Whole Aircraft Drag Polar')
xlabel('Coefficient of Lift')
ylabel('Coefficient of Drag')
title('Tempest Drag Polar')
%% Boeing Whole Aircraft Drag Polar
%Coefficient of Skin Friction
BCfe = 0.0030;
%Reference Area
BSref = 511;  
%Wetted Area
BSwet = 2896.07;    
%Oswald's Efficiency Equation
Be0 = 1.78*(1-0.045*B_AR^0.68)-0.64; 
%K value for Boeing
Bk1 = 1/(pi*Be0*B_AR);
%Minimum Coefficient of Drag
BCDmin = BCfe*(BSwet/BSref);
%Going through the 2D coefficient of drag data to find the closest min
%value
[~,i1] = min(abs(BCDmin-B_Cd));
%Finding the corresponding CLmin
BCLmin = B_CL(i1); 

Bk2 = -2 * Bk1 * BCLmin; %Check later

BCD0 = BCDmin + (Bk1 * BCLmin^2);
BCD = BCD0 + (Bk1 * B_CL.^2);

%Plotting Graphs:
%Graphs: Whole Aircraft Drag Polar (True and Model) along with Wing Drag
%Polar
figure
hold on
grid on
plot(B_CL,B_CD,'LineWidth',2)
plot(BCL_Raw,BCD_Raw,'LineWidth',2)
plot(B_CL,BCD,'LineWidth',2)
hold off
legend('Approximated Wing Drag Polar',' True Whole Aircraft Drag Polar','Approximated Whole Aircraft Drag Polar')
xlabel('Coefficient of Lift')
ylabel('Coefficient of Drag')
title('Boeing Drag Polar')

%% Tempest Performance
%Tempest Weight
T_Weight = 62.78;
%Flying Altitude
T_alt = 1500;
%Calculating L/D max for Each Coefficient of Lift and Drag (3D)
T_LD = T_CL./TCD;
[T_MaxLD,Tq] = max(T_LD);
T_RawLD = TCL_Raw./TCD_Raw;
TCD0_Raw = TCD_Raw(:) - (TCL_Raw(:).^2*Tk1);
TCD0_Raw = mean(TCD0_Raw(1:14));
[T_RawMaxLD,Ti] = max(T_RawLD);

%Density at that Altitude
T_Roe =1.0581;


%Max Glide Distance for Tempest
T_max_glide = T_alt*(T_MaxLD);
RAW_T_max_glide = T_alt*(T_RawMaxLD);

%Corresponding Angle of Attack for Max Glide Distance
TMG_AoA = TAoA(Tq); 

%Cl = L/q*s q = 1/2pv^2
% v = sqrt(2*L/(Cl*roe*s))

%Calculating Glide Velocity
TMG_V = sqrt((2*T_Weight/(T_CL(Tq)*T_Roe*TSref)));
TMG_VTrue = sqrt((2*T_Weight/(TCL_Raw(Ti)*T_Roe*TSref)));


%Max Power Range (Same as glide)
TPR_AoA = TAoA(Tq);
TPR_V = sqrt((2*T_Weight/(T_CL(Tq)*T_Roe*TSref)));
TPR_VTrue = sqrt((2*T_Weight/(TCL_Raw(Ti)*T_Roe*TSref)));

%Max Power Endurance
% CL = sqrt((cd0*3)/k1) 
TPE_CL = sqrt((TCD0*3)/Tk1);
TPE_CLTrue = sqrt((TCD0_Raw*3)/Tk1);

% v = sqrt((2*L)/CL*roe*S)
TPE_V = sqrt((2*T_Weight/(TPE_CL*T_Roe*TSref)));
TPE_VTrue = sqrt((2*T_Weight/(TCL_Raw(Ti)*T_Roe*TSref)));

%Find indicies for 
[~,TPE_i] = min(abs(TPE_CL - T_CL));
[~,TPE_iTrue] = min(abs(TPE_CLTrue - T_CL));
%Angle of attack
TPE_AoA = TAoA(TPE_i);
TPE_AoATrue = TAoA(TPE_iTrue);



%% Boeing Performance
%Boeing Weight
B_Weight = 3.559e6;
%Flying Altitude
B_alt = 10500;

%Calculating L/D max for Each Coefficient of Lift and Drag (3D)
B_LD = B_CL./BCD;
BCD0_Raw = BCD_Raw(:) - (BCL_Raw(:).^2*Bk1);
BCD0_Raw = mean(BCD0_Raw(1:19));
[MaxB_LD,Bq] = max(B_LD);
B_RawLD = BCL_Raw./BCD_Raw;
[B_RawMaxLD,Bi] = max(B_RawLD);

%Density at that Altitude
B_Roe =0.38857;

%Max Glide for the Boeing
B_max_glide = B_alt*(MaxB_LD);
RAW_B_max_glide = B_alt*(B_RawMaxLD);

%Cl = L/q*s q = 1/2pv^2
% v = sqrt(2*L/(Cl*roe*s))

%Corresponding Angle of Attack for Max Glide Distance
BMG_AoA = BAoA(Bq);
BMG_V = sqrt((2*B_Weight/(B_CL(Bq)*B_Roe*BSref)));
BMG_VTrue = sqrt((2*B_Weight/(BCL_Raw(Bi)*B_Roe*BSref)));

%Max Power Range
%Cl = sqrt((CD0*1/3)/k1)
BPR_CL = sqrt((BCD0*(1/3))/Bk1);
BPR_CLTrue = sqrt((BCD0_Raw*(1/3))/Bk1);

% v = sqrt((2*L)/CL*roe*S)
BPR_V = sqrt((2*B_Weight/(BPR_CL*B_Roe*BSref)));
BPR_VTrue = sqrt((2*B_Weight/(BCL_Raw(Bi)*B_Roe*BSref)));

%Max Power Endurance (Same as glide)
BPE_AoA = BAoA(Bq);
BPE_V = sqrt((2*B_Weight/(B_CL(Bq)*B_Roe*BSref)));
BPE_VTrue = sqrt((2*B_Weight/(BCL_Raw(Bi)*B_Roe*BSref)));

