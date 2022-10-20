%% Read In and Assigning:
%Read in the Data
table_Tempest1 = readtable("TempestData.xlsx");
table_Tempest2 = readtable("TempestData2.xlsx");
table_Boeing1 = readtable("BoeingData.xlsx");
table_Boeing2 = readtable("BoeingData2.xlsx");


%% TEMPEST CALCULATIONS: SECTION 1

%Angle of Attack Data 2D:
Tempest_1a = table_Tempest1{:,1};

%Coefficient of Lift 2D:
Tempest_1cl = table_Tempest1{:,2};

%Coefficient of Drag 2D:
Tempest_1cd = table_Tempest1{:,3};
% Calculate 3D lift slope
%Calculating a0
plot(Tempest_1a,Tempest_1cl)
[m,~] = polyfit(Tempest_1a(1:14),Tempest_1cl(1:14),1);
%xline(-2);
%xline(8);
%Calculating m3D
a0_T = m(1);
AR_T = 16.5;
e_T = 0.9;
m_3D = (a0_T)/(1+(57.3*a0_T)/(pi*e_T*AR_T));

%Obtaining the x intercept:
x_T = (-1)*m(2)/m(1);

%Obtaining CL:
CL_T = m_3D*(Tempest_1a(1:14)-x_T);

% Plotting 2D AND 3D Lift Coefficient Graphs: Tempest
plot(Tempest_1a(1:end),Tempest_1cl(1:end));
title("Coefficient of Lift against Angle of Attack for Tempest");
xlabel("Angle of Attack (degrees)");
ylabel("Coefficient of Lift");
hold on 
%origin on
grid on
plot(Tempest_1a(1:14),CL_T,'LineWidth',2);
legend("2D Lift Data","3D Lift Data");
%% BOEING CALCULATIONS: SECTION 1

%Angle of Attack 2D:
Boeing_1a = table_Boeing1{:,1};

%Coefficient of Lift 2D:
Boeing_1cl = table_Boeing1{:,2};

%Coefficient of Drag 2D:
Boeing_1cd = table_Boeing1{:,3};


% Calculate 3D lift slope
%Calculating a0
%plot(Boeing_1a,Boeing_1cl)
[m,~] = polyfit(Boeing_1a(1:17),Boeing_1cl(1:17),1);
%xline(-5)
%xline(11)

%Calculating m3D
a0_B = m(1);
AR_B = 7;
e_B = 0.9;
m_3D2 = (a0_B)/(1+(57.3*a0_B)/(pi*e_B*AR_B));

%Obtaining the x intercept:
x_B = (-1)*m(2)/m(1);

%Obtaining CL:
CL_B = m_3D2*(Boeing_1a(1:17)-x_B);

% Plotting Final Graphs:
hold off
figure
plot(Boeing_1a(1:end),Boeing_1cl(1:end));
hold on
plot(Boeing_1a(1:17),CL_B);
title("Coefficient of Lift against Angle of Attack for Boeing");
xlabel("Angle of Attack (degrees)");
ylabel("Coefficient of Lift")
legend("2D","3D");
hold off
%% TEMPEST CALCULATIONS: SECTION 2 (Indices are 1:14)
%Plotting true data  (CL on x axis and CD on y axis)
figure
plot(table_Tempest2{1:end,2},table_Tempest2{1:end,3})

%Plotting 3D approximations

%Wing Drag Polar
kT_1 = 1/(pi*e_T*AR_T);
CD_WDT = table_Tempest1{1:14,3} + (CL_T).^2*kT_1;

%Whole Aircraft Drag Polar
% Oswald's Efficiency Formula
% Tempest
e0_T = 1.78*(1-0.045*(AR_T)^(0.68))-0.64;
%CD_ADT = Parasite Drag + (CL_T)^2*(1/(pi*AR_T*e0_T))
%% SECTION 2: Boeing (Indices are 1:17)
%Plotting true data (CL on x axis and CD on y axis)
figure
plot(table_Boeing2{1:end,1},table_Boeing2{1:end,2});

%Plotting 3D approximations

%Wing Drag Polar
kB_1 = 1/(pi*e_B*AR_B);
CD_WDB = table_Boeing1{1:17,3}+(CL_B).^2*kB_1;


%Whole Aircraft Drag Polar
% Oswald's Efficiency Formula
% Tempest
e0_B = 1.78*(1-0.045*(AR_B)^(0.68))-0.64;
%CB_ADT = Parasite Drag + (CB_T)^2*(1/(pi*AR_B*e0_B))


