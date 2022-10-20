 %TBM Cutting Head Torque model 2022
clc;
close all;
clear;
tic
load SoilMonteCarlo

%% Geometric properties of the cutting head

cuttingHead.d = 0.58;                               %diameter of the cutting head [m]
cuttingHead.OR = 0.3;                               %opening ratio of the cutting head
cuttingHead.Midpoint = sqrt((cuttingHead.d/2)^2/2); %Applied location of the friction force on the main head
cuttingHead.A = (cuttingHead.d/2)^2*pi*(1-cuttingHead.OR); %Area of the front plate of the cutting head
cuttingHead.SpikeRad = 0.04;                        %Radius of the center spike [m]
cuttingHead.SpikeHeight = 0.06;                     % height of the spike
cuttingHead.SpikeAngle = atan(cuttingHead.SpikeRad/cuttingHead.SpikeHeight); %Angle of the spike relative to normal
cuttingHead.SpikeArea = cuttingHead.SpikeRad^2*pi;
cuttingHead.SpikeSA = 0.00710861;                   %surface area of the spike
cuttingHead.TeethArea = 0.00118064;                 %Base Area of each cutting tooth [m]
cuttingHead.InnerTeethCSArea = 0.00054229;          %Cross Sectional Area of the inner teeth
cuttingHead.TeethPenetration = 0.03;                %Penetration depth of the cutting teeth [m]
cuttingHead.TeethSA = 0.0035591;                    % Exposed Surface Area of each cutting tooth
cuttingHead.TeethAngle = 45;                        %cutting angle of each tooth
cuttingHead.TeethDist = [0.082,0.145,0.208];        %Radius of each row of teeth
cuttingHead.TeethLocationNum = [4,6,8];             %Amount of teeth in each row
cuttingHead.InnerTeethNum = sum(cuttingHead.TeethLocationNum); %Number of inner teeth on the cutting head
cuttingHead.OuterTeethNum = 12;                     %Number of outer angled teeth
cuttingHead.OuterTeethAngle = 35;                   %angle of the outer cutting teeth [degrees]
cuttingHead.OuterTeethSA = 0.00564661;              %Exposed Surface area of outer cutting teeth [m^2]
cuttingHead.OuterTeethCSArea = cuttingHead.InnerTeethCSArea + 3.66903*10^(-4); %Exposed Cross Sectional area of the outer teeth
cuttingHead.A_HeadBare = cuttingHead.A - cuttingHead.TeethArea * (cuttingHead.InnerTeethNum+cuttingHead.OuterTeethNum) - cuttingHead.SpikeArea; % Area of flat cutting head
cuttingHead.RPM = 2500/97;                          %RPM of the cutting head 2500 RPM/ 1:97 gear ratio
cuttingHead.Thickness = 0.0254;                     %Plate thickness of the cutting head [m] 
cuttingHead.Volume = 884;                           %Volume of the cutting head in m^3
cuttingHead.Density = 7833.413;                     %Density of Cutting Head Material kg/m^3

TorqueMotor = 77;                                   %Nm
TorqueOut = TorqueMotor * 71;                       %Nm


%% Geometric Properties of new cutting head
cuttingHead2.d = 0.58;                              %diameter of the cutting head [m]
cuttingHead2.OR = 0.335;                            %opening ratio of the cutting head
cuttingHead2.Midpoint = sqrt((cuttingHead2.d/2)^2/2); %Applied location of the friction force on the main head
cuttingHead2.A = (cuttingHead2.d/2)^2*pi*(1-cuttingHead2.OR); %Area of the front plate of the cutting head
cuttingHead2.SpikeRad = 0.045;                      %Radius of the center spike [m]
cuttingHead2.SpikeHeight = 0.06;                    % height of the spike
cuttingHead2.SpikeAngle = atan(cuttingHead2.SpikeRad/cuttingHead2.SpikeHeight); %Angle of the spike relative to normal
cuttingHead2.SpikeArea = cuttingHead.SpikeRad^2*pi;
cuttingHead2.SpikeSA = 0.01060288;                  %surface area of the spike
cuttingHead2.TeethArea = 0.00119796;                %Base Area of each cutting tooth [m]
cuttingHead2.InnerTeethCSArea = 0.00066929;         %Cross Sectional Area of the inner teeth
cuttingHead2.TeethPenetration = 0.035;              %Penetration depth of the cutting teeth [m]
cuttingHead2.TeethSA = 0.0043035;                   % Exposed Surface Area of each cutting tooth
cuttingHead2.TeethAngle = 45;                       %cutting angle of each tooth
cuttingHead2.TeethDist = [0.0772,0.11,0.1533,0.1598,0.1942,0.228,0.232]; %Radius of each row of teeth
cuttingHead2.TeethLocationNum = [4,4,4,4,4,4,4];    %Amount of teeth in each row
cuttingHead2.InnerTeethNum = sum(cuttingHead.TeethLocationNum); %Number of inner teeth on the cutting head
cuttingHead2.OuterTeethNum = 12;                    %Number of outer angled teeth
cuttingHead2.OuterTeethAngle = 35;                  %angle of the outer cutting teeth [degrees]
cuttingHead2.OuterTeethSA = 0.00564661;             %Exposed Surface area of outer cutting teeth [m^2]
cuttingHead2.OuterTeethCSArea = cuttingHead.InnerTeethCSArea + 3.66903*10^(-4); %Exposed Cross Sectional area of the outer teeth
cuttingHead2.A_HeadBare = cuttingHead.A - cuttingHead.TeethArea * (cuttingHead.InnerTeethNum+cuttingHead.OuterTeethNum) - cuttingHead.SpikeArea; % Area of flat cutting head
cuttingHead2.RPM = 24;                              %RPM of the cutting head
cuttingHead2.Thickness = 0.0254;                    %Plate thickness of the cutting head [m] 
cuttingHead2.Volume = 0.979;                        %Volume of the cutting head m^3
cuttingHead2.Density = 7833.413;                    %Density of Cutting Head Material kg/m^3

%% Geometric Properties of the PDB Cutting Head
cuttingHead3.d = 0.58;                              %diameter ofa the cutting head [m]
cuttingHead3.OR = 0.42926;                          %opening ratio of the cutting head
cuttingHead3.Midpoint = sqrt((cuttingHead.d/2)^2/2);%Applied location of the friction force on the main head
cuttingHead3.A = (cuttingHead3.d/2)^2*pi*(1-cuttingHead3.OR); %Area of the front plate of the cutting head with the teeth
cuttingHead3.SpikeRad = 0.0425;                     %Radius of the center spike [m]
cuttingHead3.SpikeHeight = 0.06;                    % height of the spike
cuttingHead3.SpikeAngle = atan(cuttingHead2.SpikeRad/cuttingHead2.SpikeHeight); %Angle of the spike relative to normal
cuttingHead3.SpikeArea = cuttingHead.SpikeRad^2*pi;
cuttingHead3.SpikeSA = 0.00981718;                  %surface area of the spike
cuttingHead3.TeethArea = 0.00118066;                %Base Area of each cutting tooth [m]
cuttingHead3.InnerTeethCSArea = 0.00066929;         %Cross Sectional Area of the inner teeth
cuttingHead3.TeethPenetration = 0.035;              %Penetration depth of the cutting teeth [m]
cuttingHead3.TeethSA = 0.0043035;                   % Exposed Surface Area of each cutting tooth
cuttingHead3.TeethAngle = 45;                       %cutting angle of each tooth
cuttingHead3.TeethDist = [0.07253,0.1123,0.1475,0.1625,0.1976,0.2325,0.2275]; %Radius of each row of teeth
cuttingHead3.TeethLocationNum = [4,4,4,4,4,4,4];    %Amount of teeth in each row
cuttingHead3.InnerTeethNum = sum(cuttingHead.TeethLocationNum); %Number of inner teeth on the cutting head
cuttingHead3.OuterTeethNum = 12;                    %Number of outer angled teeth
cuttingHead3.OuterTeethAngle = 35+cuttingHead3.TeethAngle; %angle of the outer cutting teeth  tips[degrees]
cuttingHead3.OuterTeethSA = 0.00564661;             %Exposed Surface area of outer cutting teeth [m^2]
cuttingHead3.OuterTeethCSArea = cuttingHead.InnerTeethCSArea + 3.66903*10^(-4); %Exposed Cross Sectional area of the outer teeth
cuttingHead3.A_HeadBare = cuttingHead.A - cuttingHead.TeethArea * (cuttingHead.InnerTeethNum+cuttingHead.OuterTeethNum) - cuttingHead.SpikeArea; % Area of flat cutting head
cuttingHead3.RPM = 24;                              %RPM of the cutting head
cuttingHead3.Thickness = 0.0254;                    %Plate thickness of the cutting head [m] 
cuttingHead3.Volume = 0.0061074;                    %Volume of the cutting head m^3
cuttingHead3.Density = 7833.413;                    %Density of Cutting Head Material kg/m^3

%% Soil Parameters and design Depths

FS = 1.5;

%Properties of a slitly clay mixture assuming an internal friction angle of %25 degrees and density of 1760 kg/m^3

%WE ARE NOW FAT CLAY (High Plasticity) or LEAN CLAY Cleyey Sand

%From Lean Clay data results and the Geotechnical Report Provided

%Lean Clay https://www.researchgate.net/figure/Basic-properties-of-the-lean-clay-soil_tbl2_339499761

% Calculating the equivelant fluid density of the soil from the TBC Boring
% logs, 6-8 ft depth
BoringLog.dryWeight = [96.2;101.1;100.2]*157.0874606377; %Dry soil unit weight at 6-8 ft from page 5 of report N/m^3
BoringLog.WaterContent = [0.153;0.177;0.167]; %Water content percentage of the soil
BoringLog.WetWeightAll = BoringLog.dryWeight.*(1+BoringLog.WaterContent); %Total weight of the soil in the boring logs
BoringLog.gamma = mean(BoringLog.WetWeightAll); %Average soil unit weight of the soil N/m^3
BoringLog.gammaStd = std(BoringLog.WetWeightAll); %Standard deviation of the unit weight of the soil N/m^3
BoringLog.c = 3447.38; %Cohesion constant from the boring logs for 6-8 ft depth page 5 of report Pa
BoringLog.mu = 0.55; %Friction coefficient between cohesive soils and steel
BoringLog.phi = 36.5; %Friction/slop angle of soil in degrees


softClay.gamma =   16.603*1000; %For Lean clay   % Soil unit weight (approx. for sand/silt/clay from mathalino.com) [N/m3]
softClay.mu = .55;              % max dynamic friction coeff for sand/silt/clay combo
softClay.c = 16547;  %From Boring Log %75000;             %cohesion constant of Inorganic clays, silty clays, sandy clays  - compacted, From geotechdata, [Pa]
softClay.phi = 36.5;              %Friction/slope angle
softClay.Ka = (1-sind(softClay.phi))/(1+sind(softClay.phi)); %Lateral pressure coefficient. Using active value, via robbinsfoundationsystems, using soil friction angle of 30 degrees
% 
% caliche.sigma = 11.9*10^6; %ultimate strength of caliche [Pa]
% caliche.c = 0; %cohesionless soil
% caliche.mu = 0.3;
% 
% mixture.rho_g = 100000;
% mixture.mu = 0.59;
% mixture.c = 105000;
% 
% %Sand assuming a density of 1555 kg.m^3 and an internal friction angle of
% %40 degrees
% sand.rho_g = 70154.3;
% sand.mu = 0.54;
% sand.c = 0; %Sand is a cohesionless soil no cohesion parameter

%Tunnel Properties
TunnelLength = 30;              %Length of the tunnel to be dug in m
MaxDepth = 1.5;                 %max depth of the top of the TBM m
cuttingHead.DesignDepth = MaxDepth + cuttingHead.d/2 + 0.01; %depth of the center of the TBM at the max depth

%% Calculating Soil properties
[TBCSoilProperties,TBCCoefficients] = ModelSoil(BoringLog,cuttingHead3); %Creating a model of the soil for the TBC Boring logs
[ClaySoiProperties,ClayCoefficients] = ModelSoil(softClay,cuttingHead3);

%% Solving for torque and penetration force for the Boring Log Soil on the cutting head at design depth of 1.8 m

F_ThrustTBC = calcNetThrust(cuttingHead3,TBCSoilProperties,BoringLog,0.85);
TBCTorque = CalcCuttingTorque(cuttingHead3,BoringLog,TBCSoilProperties,0.85,FS);

TBCSpeed = cuttingHead3.RPM*1*cuttingHead3.TeethPenetration;
TBCCompletionTime = TunnelLength/TBCSpeed;
format shortG
fprintf("****************For the TBC Site via Boring Logs****************\n")
fprintf("Maximum Depth: %1.2f m\n",cuttingHead.DesignDepth)
fprintf("Passive Lateral Soil Pressure at max depth: %.2f Pa\n",TBCSoilProperties.sigmaH_p(181))
fprintf("Vertical Soil Pressure at max depth: %.2f Pa\n",TBCSoilProperties.sigmaV(181))
fprintf("Forward Net Thrust Required at %1.2f m Depth: %.2f N\n",cuttingHead.DesignDepth,F_ThrustTBC(181));
fprintf("Torque Required at %1.1f m Depth: %.2f Nm\n",cuttingHead.DesignDepth,TBCTorque(181));
fprintf("Dig Speed: %.2f m/min\n",TBCSpeed);
fprintf("Time to Completion: %.2f min\n\n\n",TBCCompletionTime);

%% Solving for torque and penetration force for silty clay on the cutting head at design depth of 1.8 m
% 
% F_ThrustClay = calcThrust(cuttingHead3,ClaySoiProperties,softClay,0.75);
% ClayTorque = CalcTorque(cuttingHead3,softClay,ClaySoiProperties,0.75,FS);
% 
% ClaySpeed = cuttingHead3.RPM*1*cuttingHead3.TeethPenetration;
% ClayCompletionTime = TunnelLength/ClaySpeed;
% format shortG
% fprintf("Passive Lateral Soil Pressure at design depth: %f\n",ClaySoiProperties.sigmaH_p(181))
% fprintf("Requirements to fully penetrate silty clay:\n");
% fprintf("Forward Net Thrust Required at %f m Depth: %f N\n",cuttingHead.DesignDepth,F_ThrustClay.Total(181));
% fprintf("Torque Required at %f m Depth: %f Nm\n",cuttingHead.DesignDepth,ClayTorque.Design(181));
% fprintf("This is equal to %f Orange Subaru Crosstreks\n",ClayTorque.Design(181)/200);
% fprintf("Dig Speed: %f m/min\n",ClaySpeed);
% fprintf("Time to Completion: %f min\n",ClayCompletionTime);

% ClayRecessedSavings = (1-ClayTorque.Total(181)/ClayTorque.NonRecessed(181))*100;
% fprintf("Torque Savings from Recessed Cutting Head: %f %%\n",ClayRecessedSavings);


%% Creating graph modeling necessary torque force vs depth for silty clay


figure()
hold on
plot(TBCSoilProperties.depths,TBCSoilProperties.sigmaV,"LineWidth",1.5)
plot(TBCSoilProperties.depths,TBCSoilProperties.sigmaH_p,"LineWidth",1.5)
plot(TBCSoilProperties.depths,TBCSoilProperties.sigmaH_a,"LineWidth",1.5)
plot(TBCSoilProperties.depths,TBCSoilProperties.sigmaH_o,"LineWidth",1.5)
%plot(TBCSoilProperties.depths,TBCSoilProperties.sigmaPlate,"LineWidth",1.5)
%plot(TBCSoilProperties.depths,TBCSoilProperties.sigmaInnerTooth,"LineWidth",1.5)
%plot(TBCSoilProperties.depths,TBCSoilProperties.sigmaOuterTooth,"LineWidth",1.5)
xline(cuttingHead.DesignDepth,"k--","LineWidth",1.5)
title("Strength of Soils along different planes for the TBC Boring Log","FontSize",20)
xlabel("Depth [m]","FontSize",18)
ylabel("Soil Strength [N/m^2]","FontSize",18)
%yticklabels(0:50000:300000)
%,'Strength Normal to Plate','Inner Teeth Shear Strength','Outer teeth Shear Strength',
legend('Vertical Soil Pressure','Passive Lateral Soil Pressure','Active Lateral Soil Pressure','At Rest Lateral Soil Pressure',"Max Expected Depth","Location","NW","FontSize",16)
grid on
grid minor
ax = gca;
ax.FontSize = 16; 
hold off
% figure()
% hold on
% plot(ClaySoiProperties.depths,ClayTorque.InnerTeeth,"k--","LineWidth",1.5)
% plot(ClaySoiProperties.depths,ClayTorque.OuterTeeth,"b--","LineWidth",1.5)
% plot(ClaySoiProperties.depths,ClayTorque.Head,"Color",'#77AC30',"LineStyle","--","LineWidth",1.5)
% plot(ClaySoiProperties.depths,ClayTorque.Design/1.5,"r--","LineWidth",1.5)
% plot(ClaySoiProperties.depths,ClayTorque.Design,"Color",'#A2142F',"LineWidth",1.5)
% xline(cuttingHead.DesignDepth,"k--","LineWidth",1.5)
% title("Contributors to Overall Torque","FontSize",20)
% xlabel("Depth [m]","FontSize",18)
% ylabel("Torque [Nm]","FontSize",18)
% grid on
% legend('Inner Teeth','outer Teeth','Front Plate',"Total Torque","Total Torque with 50% increase","Max Expected Depth","Location","NW","FontSize",14)
% ax = gca;
% ax.FontSize = 16; 
% hold off

%% Finding Torque for diffferent penetrations

Penetrations = (0.1:0.01:1)';   % If we dont have full penetration into the soil (which we wont)
F_ThrustClayVar = zeros(length(TBCSoilProperties.depths),length(Penetrations));
F_ThrustClayVarDesignDepth  = zeros(length(Penetrations),1);
F_ThrustClayVarHalfDesignDepth = zeros(length(Penetrations),1);
F_ThrustClayVarThreeQuarterDesignDepth = zeros(length(Penetrations),1);
for i = 1:length(Penetrations)
    F_ThrustClayVar(:,i) = calcNetThrust(cuttingHead3,ClaySoiProperties,softClay,Penetrations(i));
    F_ThrustClayVarDesignDepth(i) = F_ThrustClayVar(181,i);
    F_ThrustClayVarHalfDesignDepth(i) = F_ThrustClayVar(91,i);
    F_ThrustClayVarThreeQuarterDesignDepth(i) = F_ThrustClayVar(136,i);
end

figure()
hold on
grid on
grid minor
plot(Penetrations,F_ThrustClayVarDesignDepth,"LineWidth",1.5)
plot(Penetrations,F_ThrustClayVarHalfDesignDepth,"LineWidth",1.5)
plot(Penetrations,F_ThrustClayVarThreeQuarterDesignDepth,"LineWidth",1.5)
title("Thrust Required vs penetrations","FontSize",18)
xlabel("Penetration","FontSize",16)
ylabel("Thrust Required [N]","FontSize",16)
legend("1.5 m","0.6 m","1.05 m","Location","NW","FontSize",16)
hold off


%% Plotting torque for each penetration and time to completetion

DesignDepthTorques = zeros(length(Penetrations),1);
ClayTorqueVar = zeros(length(TBCSoilProperties.depths),length(Penetrations));
HalfDesignDepthTorques = zeros(length(Penetrations),1);
ThreeQuartersDesignDepthTorques = zeros(length(Penetrations),1);
for i = 1:length(Penetrations)
    ClayTorqueVar(:,i) = CalcCuttingTorque(cuttingHead3,softClay,ClaySoiProperties,Penetrations(i),FS);
    DesignDepthTorques(i) = ClayTorqueVar(181,i);
    HalfDesignDepthTorques(i) = ClayTorqueVar(91,i);
    ThreeQuartersDesignDepthTorques(i) = ClayTorqueVar(136,i);
end

PenetrationSpeedVar = cuttingHead3.RPM*Penetrations*cuttingHead3.TeethPenetration;
CompletionTimeVar = TunnelLength./PenetrationSpeedVar;

%% Plotting Torque vs penetration


figure()

hold on
grid on

plot(Penetrations,DesignDepthTorques,"LineWidth",1.5)
plot(Penetrations,HalfDesignDepthTorques,"LineWidth",1.5)
plot(Penetrations,ThreeQuartersDesignDepthTorques,"LineWidth",1.5)
ylabel("Torque Required [N]","FontSize",16)
yyaxis right
ylim([0,120])

plot(Penetrations,CompletionTimeVar,"LineWidth",1.5)
ylabel("Time to completetion [min]","FontSize",16)
title("Torque Required vs penetrations","FontSize",16)
xlabel("Penetration","FontSize",16)
xlim([0.1,1])
legend("Torque Required 1.5 m","Torque Required 0.6 m","Torque Required 1 m","Time to Completion","Location","NW","FontSize",16)
ax = gca;
ax.FontSize = 16;
hold off

%% Performing Multidemensional montecarlo Analysis

% VarGamma = normrnd(20000,2000,[10000,1]);
% VarC = normrnd(75000,5000,[10000,1]); %Varying the cohesion constant of the soils [Pa]
% VarMu = normrnd(0.55,0.05,[10000,1]);
% VarGamma(find(VarGamma < 10000,10)) = 10000;
% 
% VarSoil.gamma = 20000; % Soil unit weight (approx. for sand/silt/clay from mathalino.com) [N/m3]
% VarSoil.mu = .57; % max dynamic friction coeff for sand/silt/clay combo
% VarSoil.c = 86000; %cohesion constant of Inorganic clays, silty clays, sandy clays of low plasticity - compacted, From geotechdata, [Pa]
% VarSoil.phi = 30; %Friction/slope angle
% VarSoil.Ka = (1-sind(softClay.phi))/(1+sind(softClay.phi)); %Lateral pressure coefficient. Using active value, via robbinsfoundationsystems, using soil friction angle of 30 degrees
% VarSoilProperties = struct('sigmaPlate', cell(length(VarGamma), 1), 'sigmaInnerTooth', cell(length(VarGamma), 1),'sigmaOuterTooth',cell(length(VarGamma), 1),'sigmaH',cell(length(VarGamma), 1),'sigmaV',cell(length(VarGamma), 1),'sigmaSpike',cell(length(VarGamma), 1),'depths',cell(length(VarGamma), 1));
% VarSoilData = struct('gamma',cell(length(VarGamma), 1),'mu',cell(length(VarGamma), 1),'c',cell(length(VarGamma), 1),'phi',cell(length(VarGamma), 1),'Ka',cell(length(VarGamma), 1));
% for i = 1:length(VarGamma)
%     VarSoil.gamma = VarGamma(i);
%     VarSoil.c = VarC(i);
%     VarSoil.mu = VarMu(i);
%     VarSoilProperties(i) =  ModelSoil(VarSoil,cuttingHead3);
%     VarSoilData(i) = VarSoil;
% end
% % multidimensional thrust data
% VarThrustData = struct('Total',cell(length(VarGamma), length(Penetrations)));
% VarThrustDataDesignDepthWhole = zeros(length(VarSoilProperties),1);
% 
% for i = 1:length(VarSoilProperties)
%     parfor j = 1:length(Penetrations)
%         VarThrustData(i,j) = calcThrust(cuttingHead3,VarSoilProperties(i),VarSoilData(i),Penetrations(j));
%         VarTorqueData(i,j) = CalcTorque(cuttingHead3,VarSoilData(i),VarSoilProperties(i),Penetrations(j),FS);
%     end
% end

clear i j l k

%% Finding Relevant Data of MonteCarlo MultiDimensional

% VarThrustMaxDepth = zeros(length(VarThrustData),91);
% VarTorqueMaxDepth = zeros(length(VarThrustData),91);
% 
% for i = 1:length(VarThrustData)
%     for j = 1:91
%         VarThrustMaxDepth(i,j) = VarThrustData(i,j).Total(181);
%         VarTorqueMaxDepth(i,j) = VarTorqueData(i,j).Design(181);
%     end
% end
%save("SoilMonteCarlo.mat","VarTorqueMaxDepth","VarThrustMaxDepth")

VarTorqueDataDesignDesign = zeros(length(VarThrustMaxDepth),1);
VarThrustDataDesignDepthDesign = zeros(length(VarThrustMaxDepth),1);
for i = 1:length(VarThrustMaxDepth)
    VarThrustDataDesignDepthDesign(i) = VarThrustMaxDepth(i,81);
    VarTorqueDataDesignDesign(i) = VarTorqueMaxDepth(i,81);
end

clear i j

%% Mean and Standard deviation of histogram montecarlo

VarTorqueMean = mean(VarTorqueDataDesignDesign);
VarTorqueStd = std(VarTorqueDataDesignDesign);
VarThrustMean = mean(VarThrustDataDesignDepthDesign);
VarThrustStd = std(VarThrustDataDesignDepthDesign);

%% Plotting Histogram of multidimensional montecarlo data

figure()
hold on
histogram(VarThrustDataDesignDepthDesign,200)
xlabel("Thrust N","FontSize",18)
ylabel("Occurances","FontSize",18)
title("Monte Carlo Simulation of Thrust","FontSize",18)
xline(VarThrustMean-VarThrustStd,"k--","LineWidth",2,"Label","Uncertainty","FontSize",16,"LabelHorizontalAlignment","center")
xline(VarThrustMean+VarThrustStd,"k--","LineWidth",2,"Label","Uncertainty","FontSize",16,"LabelHorizontalAlignment","left")
xline(VarThrustMean,"r--","LineWidth",2,"Label","Mean","FontSize",16)
ax = gca;
ax.FontSize = 16;

% figure()
% histogram(VarThrustDataDesignDepthWhole,100)
% xlabel("Thrust")

figure()
hold on
histogram(VarTorqueDataDesignDesign,200)
xlabel("Torque Nm","FontSize",18)
ylabel("Occurances","FontSize",18)
title("Monte Carlo Simulation of Torque","FontSize",18)
xline(VarTorqueMean-VarTorqueStd,"k--","LineWidth",2,"Label","Uncertainty","FontSize",16)
xline(VarTorqueMean+VarTorqueStd,"k--","LineWidth",2,"Label","Uncertainty","FontSize",16)
xline(VarTorqueMean,"r--","LineWidth",2,"Label","Mean","FontSize",16)
ax = gca;
ax.FontSize = 16;
hold off

%% Calculating Thrust and torque profile for the dig.


path = readmatrix("pathPlanNew.txt");
path(:,2) = abs(path(:,2));
TorqueMax = 4000; % Max rated torque of the gearbox with our use case Nm

[pathTorques,pathThrusts,pathPenetrations,DigTime,RemovalSpeed] = modelPath(path,ClaySoiProperties,TorqueMax,ClayTorqueVar,F_ThrustClayVar,Penetrations,cuttingHead3);
figure()
subplot(4,1,1)
Title = "Path Data for a Max Torque of "+compose("%5.2f",TorqueMax)+" Nm. Dig Time of "+compose("%5.2f",DigTime)+" min.";
sgtitle(Title,"FontSize",20)
colormap jet
color_line3(path(:,1),-1*path(:,2),zeros(length(path),1),pathTorques,"LineWidth",4);
colorbar
%ylim([-5,0])
title("Torques over the Dig Path")
xlabel("Horizontal Distance","FontSize",16)
ylabel("Depth","FontSize",16)
ax = gca;
ax.FontSize = 16;
Bar = colorbar;
colorTitleHandle = get(Bar,'Title');
set(colorTitleHandle ,'String',"Torque Outputed");

subplot(4,1,2)
colormap jet
color_line3(path(:,1),-1*path(:,2),zeros(length(path),1),pathThrusts,"LineWidth",4);
colorbar
%ylim([-5,0])
title("Thrusts over the Dig Path")
xlabel("Horizontal Distance","FontSize",16)
ylabel("Depth","FontSize",16)
ax = gca;
ax.FontSize = 16;
Bar = colorbar;
colorTitleHandle = get(Bar,'Title');
set(colorTitleHandle ,'String',"Thrust Outputed");
Bar.TickLabels = (5000:5000:25000);

subplot(4,1,3)
colormap jet
color_line3(path(:,1),-1*path(:,2),zeros(length(path),1),pathPenetrations,"LineWidth",4);
colorbar
%ylim([-5,0])
title("Teeth Penetration over the Dig Path")
xlabel("Horizontal Distance","FontSize",16)
ylabel("Depth","FontSize",16)
ax = gca;
ax.FontSize = 16;
Bar = colorbar;
colorTitleHandle = get(Bar,'Title');
set(colorTitleHandle ,'String',"Teeth Penetration");

subplot(4,1,4)
colormap jet
color_line3(path(:,1),-1*path(:,2),zeros(length(path),1),[RemovalSpeed;RemovalSpeed(end)],"LineWidth",4);
colorbar
%ylim([-5,0])
title("Soil Removal Rate over the Dig Path")
xlabel("Horizontal Distance","FontSize",16)
ylabel("Depth","FontSize",16)
ax = gca;
ax.FontSize = 16;
Bar = colorbar;
colorTitleHandle = get(Bar,'Title');
set(colorTitleHandle ,'String',"Soil Removal Rate [gal/min]");

toc
%% Helper Functions


%Finding the optimal torque and thrust at each element along a path
function [PathTorque,PathThrust,PathPenetration,DigTime,RemovalSpeed,Speeds] = modelPath(path,soilProperties,TorqueLimit,TorqueData,ThrustData,penetrations,cuttingHead)
    
    %Allocate Vectors
    PathPenetration = zeros(length(path),1);
    PathTorque = zeros(length(path),1);
    PathThrust = zeros(length(path),1);
    
    %Going through the length of the path to find the torque at each point in
    %the path
    for i = 1:length(path)
        %Finding the corresponding index for the path depth at this point
        depth = path(i,2);
        DepthIdx = find(soilProperties.depths > depth,1,"first");
        ValTooBig = true;
        j = length(penetrations);
        while ValTooBig == true
            TempTorque = TorqueData(DepthIdx,j);
            if TempTorque <= TorqueLimit %If the rorque works
                PathPenetration(i) = penetrations(j);
                PathTorque(i) = TempTorque;
                PathThrust(i) = ThrustData(DepthIdx,j);
                ValTooBig = false;
                j=j+1;
            else
                %Repeat?
                j = j-1;
            end
        end
    end
    % Speeds = RPM*PathPen*Pen
    Speeds = cuttingHead.RPM*PathPenetration(2:end)*cuttingHead.TeethPenetration; % Speed of the cutting head between each path point
    Dist = (diff(path(:,1)).^2+diff(path(:,2)).^2).^(1/2);
    Times = Dist./Speeds;
    DigTime = sum(Times)*1.5;
    RemovalSpeed = Speeds*0.3^2*pi*264.172;

end






