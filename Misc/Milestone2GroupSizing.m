clc
clear

%% Data Extraction
data = importdata('app2data.mat');
AoA = data(:,1); %Angle of Attack (deg)
Cl = data(:,2); %2-D Coefficient of lift
Cd = data(:,3); %2-D Coefficient of Drag
Clmax = max(Cl); %Maximum Cl

%% Constants
g = 9.81; %Gravitational Acceleration (m/s^2)
Cfe = 0.003; %Skin Friction
e = 0.9; %Span Efficiency
p = 1.0581; %Air Density (kg/m^3)
h = 7; %Launch Height (m)
payload = 0.133*g; %Mass of camera (N)
ballast = 0.05*g; %Ballast weight (N)
Sother = 0.214; %Area of Aircraft excluding wings (m^2)
AR = 3; %Aspect Ratio
k = 1/(pi*e*AR);
e0 = 1.78*(1-0.045*AR^0.68)-0.64; %Oswald Efficiency 

%% Wing Loading
Sref = 0.0004:0.00001:5; %Planar Wing Area (m^2)
WAfoam = 0.295*g; %Wfoam/Afoam (N/m^2)
Wload = (WAfoam*(2*Sref+Sother) + payload + ballast)./Sref; %Wing load (N/m^2)

%% Varitation of Weight
Wto = Wload.*Sref; %Takeoff Weight (N)

%% Glide Range
Swet = (2*Sref) + Sother; %Aircraft Surface Area (m^2)
CD0 = Cfe*(Swet./Sref); %Parasite Drag Coefficient
CL_r = sqrt(CD0/k); %Coefficient of Lift for Glide Range
Grange = (sqrt(CD0/k)./(2*CD0))*h; %Glide Range (m)

%% Max Endurance
CD_e = 4.*CD0; %Coefficient of Drag for Endurance
Ve = sqrt(Wload./(0.5.*p.*sqrt(3.*CD0.*pi.*e0.*AR))); %Endurance Velocity (m/s)
qe = 0.5.*p.*Ve.^2; %Dynamic Pressure for Endurance (Pa)
E = Wload.*h./(CD_e.*qe.*Ve); %Endurance (s)
CL_e = sqrt(3.*CD0*pi*e0*AR); %Coefficient of Lift for Glide Endurance

%% Velocity Required for Max E and R
Vr = sqrt(Wload./(sqrt(CD0.*pi.*e0.*AR).*0.5.*p)); %Range Velocity (m/s)
Vstall = sqrt(2*Wto./(Clmax.*p.*Sref)); %Stall Velocity (m/s)

%% Wing Geometery
b = sqrt(AR.*Sref); %Wing Span (m)

%% Variation of Aircraft Cost
Cost = ((Wto-payload)/g)*1000; %Cost of Aircraft ($)

%% Subplots
subplot(2,3,1)
hold on
yyaxis left
plot(Wload,E,'LineWidth',2)
yline(7,'g--','LineWidth',2)
yline(10,'r--','LineWidth',2)
xlabel('Wing Load (N/m^2)')
ylabel('Glide Endurance (s)')

yyaxis right
plot(Wload,CL_e,'LineWidth',2)
ylabel('Coefficient of Lift')

legend('Emax Achieved','Min Emax','Max Emin','CL Req','Location','North')
title('Max Endurance & CL vs Wing Load')
xlim([10 45])
hold off

subplot(2,3,2)
hold on
yyaxis left
plot(Wload,Grange,'LineWidth',2)
yline(70,'g--','LineWidth',2)
yline(100,'r--','LineWidth',2)
xlabel('Wing Load (N/m^2)')
ylabel('Glide Range (m)')

yyaxis right
plot(Wload,CL_r,'LineWidth',2)
ylabel('Coefficient of Lift')

legend('Rmax Achieved','Min Rmax','Max Rmin','CL Req','Location','North')
title('Max Glide Range & CL vs Wing Load')
xlim([10 45])
hold off

subplot(2,3,3)
plot(Wload,Wto,'LineWidth',2)
xlabel('Wing Load (N/m^2)')
ylabel('Aircraft Takeoff Weight (N)')
title('Varitation of Weight with Wing Load')
xlim([10 45])

subplot(2,3,4)
hold on
plot(Wload,Vr,'LineWidth',2)
plot(Wload,Ve,'LineWidth',2)
plot(Wload,Vstall,'r--','LineWidth',2)
yline(7,'g--','LineWidth',2)
yline(12,'r--','LineWidth',2)
xlabel('Wing Load (N/m^2)')
ylabel('Velocity (m/s)')
legend('Max Range Velocity','Max Endurance Velocity','Stall Velocity',...
    'Minimum Velocity','Maximum Velocity')
title('Max E and Max R Velocity vs Wing Load')
xlim([10 45])
hold off

subplot(2,3,5)
hold on
yyaxis left
plot(Wload,b,'LineWidth',2)
yline(1,'g--','LineWidth',2)
xlabel('Wing Load (N/m^2)')
ylabel('Wing Span (m)')

yyaxis right
plot(Wload,Sref,'LineWidth',2)
ylabel('Wing Planar Area (m^2)')

legend('Span','Max Span','Planar Area','Location','North')
title('Variation of Wing Geometery')
xlim([10 45])
hold off

subplot(2,3,6)
plot(Wload,Cost,'LineWidth',2)
xlabel('Wing Load (N/m^2)')
ylabel('Aircraft Cost ($)')
title('Variation of Cost with Wing Load')
xlim([10 45])