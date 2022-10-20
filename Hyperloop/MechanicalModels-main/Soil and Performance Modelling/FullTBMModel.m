%% Header
% Full Dynamics and Perfomance Model of the TBM for FDP 2022
clc
close all;
clear
tic
load path3dUpdated.mat
load sampleTime4Baro.mat
path3d = [downrange,crossrange,depth]; %Putting the path into a single matrix
%% Geometric Properties of the Cutting Head
cuttingHead.d = 0.58;                              %diameter ofa the cutting head [m]
cuttingHead.OR = 0.42926;                          %opening ratio of the cutting head
cuttingHead.Midpoint = sqrt((cuttingHead.d/2)^2/2);%Applied location of the friction force on the main head
cuttingHead.A = (cuttingHead.d/2)^2*pi*(1-cuttingHead.OR); %Area of the front plate of the cutting head with the teeth
cuttingHead.SpikeRad = 0.0425;                     %Radius of the center spike [m]
cuttingHead.SpikeHeight = 0.06;                    % height of the spike
cuttingHead.SpikeAngle = atan(cuttingHead.SpikeRad/cuttingHead.SpikeHeight); %Angle of the spike relative to normal
cuttingHead.SpikeArea = cuttingHead.SpikeRad^2*pi;
cuttingHead.SpikeSA = 0.00981718;                  %surface area of the spike
cuttingHead.TeethArea = 0.00118066;                %Base Area of each cutting tooth [m]
cuttingHead.InnerTeethCSArea = 0.00066929;         %Cross Sectional Area of the inner teeth
cuttingHead.TeethPenetration = 0.035;              %Penetration depth of the cutting teeth [m]
cuttingHead.TeethSA = 0.0043035;                   % Exposed Surface Area of each cutting tooth
cuttingHead.TeethAngle = 45;                       %cutting angle of each tooth
cuttingHead.TeethDist = [0.07253,0.1123,0.1475,0.1625,0.1976,0.2325,0.2275]; %Radius of each row of teeth
cuttingHead.TeethLocationNum = [4,4,4,4,4,4,4];    %Amount of teeth in each row
cuttingHead.InnerTeethNum = sum(cuttingHead.TeethLocationNum); %Number of inner teeth on the cutting head
cuttingHead.OuterTeethNum = 12;                    %Number of outer angled teeth
cuttingHead.OuterTeethAngle = 35+cuttingHead.TeethAngle; %angle of the outer cutting teeth  tips[degrees]
cuttingHead.OuterTeethSA = 0.00564661;             %Exposed Surface area of outer cutting teeth [m^2]
cuttingHead.OuterTeethCSArea = cuttingHead.InnerTeethCSArea + 3.66903*10^(-4); %Exposed Cross Sectional area of the outer teeth
cuttingHead.A_HeadBare = cuttingHead.A - cuttingHead.TeethArea * (cuttingHead.InnerTeethNum+cuttingHead.OuterTeethNum) - cuttingHead.SpikeArea; % Area of flat cutting head
cuttingHead.RPM = 30;                              %RPM of the cutting head
cuttingHead.Thickness = 0.0254;                    %Plate thickness of the cutting head [m] 
cuttingHead.Volume = 0.0061074;                    %Volume of the cutting head m^3
cuttingHead.Density = 7833.413;                    %Density of Cutting Head Material kg/m^3

%% Geometric Properties of the TBM
TBM.mFront = 430; %Mass of the front Section of the TBM (kg)
TBM.SRShellCSAreaHorz = 0.13*0.3; %Horizontal Cross Sectional Area of the soil removal shell (m^2)
TBM.SRShellCSAreaVert = TBM.SRShellCSAreaHorz; %Vertical Cross Sectional Area of the soil removal Shell (m^2)
TBM.GripShellCSAreaHorz = 0.107814*2; %Horizontal Cross Sectional Area of the shell around the gripper and gearbox section (m^2)
TBM.GripShellCSAreaVert = 0.12727*2; %Vertical Cross Sectional Area of the shell around the gripper and gearbox section (m^2)
TBM.FrontGripCSAreaVert = 0.065975; %Vertical Cross Sectional Area of the Front Gripper (m^2)
TBM.FrontGripCSAreahorz = 0.168; %Horizontal Cross Sectional Area of the front Gripper (m^2)
TBM.LSRChamber = 0.13; %Length of the Soil Removal Chamber Meters
TBM.LFront = 0.77; %Length of the Front of the TBM m
TBM.LGripShell = 0.64; %Length of the Gripper and gearbok Section
TBM.Grippermu = 0.55; %Coefficient of friction between the gripper pads and the soil
TBM.HexapodShellExtCSArea = 0.491744*0.6; %Cross Sectional Area of the moving Hexapod Shell when extended
TBM.HexapodShellContCSArea = 0.288544 * 0.6; %Cross sectional area of half of the hexapod when contracted m^2
TBM.mHexapod = 115; %Mass of the hexapod in kg
TBM.TSCSArea = 0.6*1.1; %Cross Sectional Area of Tunnel Support Storage Shell m^2
TBM.LStorage = 1.1; %Length of TS Storage m 
TBM.LReleaseSection = 0.7; %Length of the Rlease arm
TBM.ReleaseShellCSArea = TBM.LReleaseSection * 0.56; %Cross Sectional Area of 
TBM.LBackGrip = 0.55; %Length of the back gripper housing m
TBM.LBack = 2.35; %Length of the Back of the TBM m
TBM.BackGripCSAreaVert = 0.1061; %Vertical Cross Sectional area of the back gripper shell m^2
TBM.BackGripCSAreaHorz = 0.0397; %Horizontal Cross Sectional area of the back gripper shell m^2
TBM.TarpSteelmus = 0.2; %Coefficient of Friction between Steel and Polyethelyne
TBM.mBack = 533; %Mass of the back section in kg

%% Data about the Actuators
PA13.maxSpeed = 0.51* 1.524; %No Load speed of the actuator m/min
PA13.minSpeed = 0.28 * 1.524; %Full load speed of the actuator m/min
PA13.maxForce = 60000; %Max Propulsive Force of All 6 Actuators
PA13.SpeedConstant = (PA13.minSpeed-PA13.maxSpeed)/PA13.maxForce; %Speed constant of the hexapod (m/min)/N
PA13.gripperStokeTime = 1; %Time for the grippers to extend (s)
PA13.StrokeLength = 8*0.0254; %Stroke Length of the Actuator (m)
PA13.DutyCycle = 1; %On off duty cycle of 0.5
PA13.sampleTime = 0; %Sample time required for the barometers (s)

%% key constants
FS = 1.5; %Factor of Safety to use
MaxDepth = 1.5; %max depth of the top of the TBM m
cuttingHead.MaxDepth = MaxDepth + cuttingHead.d/2 + 0.01; %depth of the center of the TBM at the max depth


%% Soil Data for the TBC Boring Logs
% Calculating the equivelant fluid density of the soil from the TBC Boring
% logs, 6-8 ft depth
BoringLog.dryWeight = [96.2;101.1;100.2]*157.0874606377; %Dry soil unit weight at 6-8 ft from page 5 of report N/m^3
BoringLog.WaterContent = [0.153;0.177;0.167]; %Water content percentage of the soil
BoringLog.WetWeightAll = BoringLog.dryWeight.*(1+BoringLog.WaterContent); %Total weight of the soil in the boring logs
BoringLog.gamma = mean(BoringLog.WetWeightAll); %Average soil unit weight of the soil N/m^3
BoringLog.gammaStd = std(BoringLog.WetWeightAll); %Standard deviation of the unit weight of the soil N/m^3
BoringLog.c = 3447.38; %Cohesion constant from the boring logs for 6-8 ft depth page 5 of report Pa
BoringLog.mu = 0.5; %Friction coefficient between cohesive soils and steel kinetic
BoringLog.mus = 0.55; %Friction coefficient between cohesive soils and steel static
BoringLog.phi = 36.5; %Friction/slop angle of soil in degrees

%% Calculating Soil properties
[TBCSoilData,TBCCoefficients] = ModelSoil(BoringLog,cuttingHead); %Creating a model of the soil for the TBC Boring logs

%Max Pressure

[MaxSv,index] = max(TBCSoilData.sigmaV);
MaxShp = max(TBCSoilData.sigmaH_p);
MaxSho = max(TBCSoilData.sigmaH_o);
MaxSha = max(TBCSoilData.sigmaH_a);



%% Calculating the Propulsion Force at Various Depths and Penetrations
FricFront = calcFrontFriction(TBM,TBCSoilData,BoringLog);
FricBack = calcBackFriction(TBM,TBCSoilData,BoringLog);

% Calculating the Propulsion force for each depth and penetration for the
% dig sim
Penetrations = (0.01:0.005:1)';   % Array of teeth Penetrations
F_ThrustVar = zeros(length(FricFront.Kinetic),length(Penetrations));
NetThrustVar = zeros(length(FricFront.Kinetic),length(Penetrations));
for i = 1:length(Penetrations)
    NetThrustTemp = calcNetThrust(cuttingHead,TBCSoilData,BoringLog,Penetrations(i));
    NetThrustVar(:,i) = NetThrustTemp(1:length(FricFront.Kinetic));
    F_ThrustVar(:,i) = NetThrustVar(:,i) + FricFront.Kinetic;
end
clear NetThrustTemp i

%% Calculating the Teeth Penetration Profile over the dig
DigPenetration = calcPenetration(cuttingHead,F_ThrustVar,PA13,Penetrations,TBCSoilData.depths,path3d); %Penetrations over the course of the dig


%% Calculating the Torques and Thrusts over the Dig
TorqueDig = zeros(length(DigPenetration),1);

path3d(:,3) = abs(path3d(:,3)); %Fixing the coordinate frame
for i = 1:length(DigPenetration)
    DepthIdx = find(TBCSoilData.depths > path3d(i,3),1,"first");
    PenIdx = find(Penetrations == DigPenetration(i),1);
    SoilDataDepth.sigmaV = TBCSoilData.sigmaV(DepthIdx);
    SoilDataDepth.sigmaH_p = TBCSoilData.sigmaH_p(DepthIdx);
    SoilDataDepth.sigmaH_o = TBCSoilData.sigmaH_o(DepthIdx);
    SoilDataDepth.sigmaInnerTooth = TBCSoilData.sigmaInnerTooth(DepthIdx);
    SoilDataDepth.sigmaOuterTooth = TBCSoilData.sigmaOuterTooth(DepthIdx);
    SoilDataDepth.sigmaSpike = TBCSoilData.sigmaSpike(DepthIdx);
    FricFrontPath.Kinetic(i) = FricFront.Kinetic(DepthIdx);
    FricFrontPath.Static(i) = FricFront.Static(DepthIdx);
    FricBackPath.Kinetic(i) = calcBackFriction(TBM,TBCSoilData,BoringLog).Kinetic(DepthIdx); %Friction Force Required to overcome to Pull the back forward
    FricBackPath.Static(i) = calcBackFriction(TBM,TBCSoilData,BoringLog).Static(DepthIdx);
    TorqueDig(i) = CalcCuttingTorque(cuttingHead,BoringLog,SoilDataDepth,DigPenetration(i),FS); %Torque at every point along the dig Nm
    ThrustDig.NetThrust(i) = NetThrustVar(DepthIdx,PenIdx);
    ThrustDig.Total(i) = ThrustDig.NetThrust(i) + FricFrontPath.Kinetic(i); %Required Thrust at each point N No FS
    ThrustDig.Design(i) = ThrustDig.Total(i)*FS; %Required Thrust at each point N with FS

end

ThrustDig.Total = ThrustDig.Total';
ThrustDig.Design = ThrustDig.Design';
ThrustDig.NetThrust = ThrustDig.NetThrust';
FricFrontPath.Kinetic = FricFrontPath.Kinetic';
FricFrontPath.Static = FricFrontPath.Static';
FricBackPath.Kinetic = FricBackPath.Kinetic';
FricBackPath.Static = FricBackPath.Static';
clear i  SoilDataDepth DepthIdx PenIdx
ThrustMax = max(ThrustDig.Total); %Maximum thrust required
ThrustMaxDesign = max(ThrustDig.Design); %Maximum thrust required
NetThrustMax = max(ThrustDig.NetThrust);
MaxDepthIdx = find(TBCSoilData.depths > max(path3d(:,3)),1,"first");
FricMax = FricFront.Kinetic(MaxDepthIdx);
FricBackMax = FricBack.Static(MaxDepthIdx);
FricBackKineticMax = FricBack.Kinetic(MaxDepthIdx);
FricFrontStaticMax = FricFront.Static(MaxDepthIdx);
MaxBackGripperForce = (ThrustMax-FricBackMax)/TBM.Grippermu/2;
MaxFrontGripperForce = (FricBackKineticMax-FricFrontStaticMax)/TBM.Grippermu/2;
TorqueMax = max(TorqueDig);
fprintf("Max Thrust Force: %5.0f N\n",ThrustMax)
fprintf("Max Thrust Force with %1.2f FS: %5.0f N\n",FS,ThrustMaxDesign)
fprintf("Max Net Thrust Force: %5.0f N\n",NetThrustMax)
fprintf("Max Frictional Force: %5.0f N\n",FricMax)
fprintf("Max Frictional Force on Back of TBM: %5.0f N\n",FricBackMax)
fprintf("Back Gripper Force Required on Single Gripper: %5.0f N \n",MaxBackGripperForce)
fprintf("Front Gripper Force Required on Single Gripper: %5.0f N \n",MaxFrontGripperForce)
fprintf("Max Torque with %1.2f FS: %4.0f Nm\n",FS,TorqueMax)

%% Calculating the time to completion of the Dig
[DigTime,VActive] = digTime(ThrustDig.Total,FricBackPath.Kinetic,path3d,PA13);
fprintf("Total Dig Time: %3.0f min \n",DigTime.TotalTime)
fprintf("Active Dig Time: %2.0f min \n",DigTime.ActiveDigTime)
fprintf("Pull Time: %2.0f min \n",DigTime.PullTime)
fprintf("Time For Grippers to Extend: %2.0f min \n",DigTime.GripperExtendTime)
fprintf("Stationary Time: %3.0f min \n",DigTime.Stationary)

%DigThrust = [struct2table(ThrustDig),array2table(FricBackPath.Kinetic,"VariableNames","Pull Thrust")];

%% Graphing the Key metrics over the dig
figure()
subplot(3,1,1)
Title = "Dynamics Model for the Dig";
sgtitle(Title,"FontSize",20)
colormap jet
color_line3(path3d(:,1),path3d(:,3),zeros(length(path3d(:,3)),1),TorqueDig,"LineWidth",4);
grid on
grid minor
colorbar
TorqueTitle = "Torques over the Dig Path";
title(TorqueTitle,"FontSize",18)
xlabel("Horizontal Distance (m)","FontSize",16)
ylabel("Depth (m)","FontSize",16)
xlim([0,path3d(end,1)+0.1])
ax = gca;
ax.FontSize = 16;
Bar = colorbar;
colorTitleHandle = get(Bar,'Title');
set(colorTitleHandle ,'String',"Torque Outputed");
ax = gca;
ax.YDir = "reverse";

subplot(3,1,2)
colormap jet
color_line3(path3d(:,1),path3d(:,3),zeros(length(path3d(:,3)),1),ThrustDig.Design,"LineWidth",4);
grid on
grid minor
colorbar
TorqueTitle = "Thrust over the Dig Path";
title(TorqueTitle,"FontSize",18)
xlabel("Horizontal Distance (m)","FontSize",16)
ylabel("Depth (m)","FontSize",16)
xlim([0,path3d(end,1)+0.1])
ax = gca;
ax.FontSize = 16;
Bar = colorbar;
colorTitleHandle = get(Bar,'Title');
set(colorTitleHandle ,'String',"Torque Outputed");
ax = gca;
ax.YDir = "reverse";

subplot(3,1,3)
colormap jet
color_line3(path3d(:,1),path3d(:,3),zeros(length(path3d(:,3)),1),DigPenetration,"LineWidth",4);
grid on
grid minor
colorbar
TorqueTitle = "Teeth Penetration over the Dig Path";
title(TorqueTitle,"FontSize",18)
xlabel("Horizontal Distance (m)","FontSize",16)
ylabel("Depth (m)","FontSize",16)
xlim([0,path3d(end,1)+0.1])
ax = gca;
ax.FontSize = 16;
Bar = colorbar;
colorTitleHandle = get(Bar,'Title');
set(colorTitleHandle ,'String',"Torque Outputed");
ax = gca;
ax.YDir = "reverse";

toc



