%Calculates the friction required to push the front section at a given
%depth
function Friction = calcFrontFriction(TBM,soilData,soilProperties)
% Inputs:
% TBM: The Geometric Properties of the TBM
% soilData: the stengths of the soil at different depths in different
% planes
% soilProperties: The parameters of the soil that are constant with depth
% FS: Desired Factor of Safety

    F_nTBM = TBM.mFront*9.81;%Normal force of the TBM
    %Getting the Necessary Soil Pressures
    P_vertUpper = [zeros(26,1);soilData.sigmaV(1:end-52)]; %Vertical Pressure on the upper half of the TBM Approximating at the centriod of a semicircle
    P_HorzUpper = [zeros(26,1);soilData.sigmaH_o(1:end-52)]; %Horizontal Pressure on upper half of TBM, using at rest pressure because soil is statonary relative to TBM
    P_HorzLower = soilData.sigmaH_o(27:end); %Horizontal Pressure on the bottom of the TBM, using at rest pressure
    P_HorzMid = soilData.sigmaH_o(1:end-26); %Horizontal Pressure at the mid point of the TBM, using at rest pressure

    %Soil Removal Shell
    F_SRVert = P_vertUpper *TBM.SRShellCSAreaVert; %Force due to vertical Pressure on soil removal chamber
    F_SRHorzU = P_HorzUpper*TBM.SRShellCSAreaHorz; %Force on upper half of SR Chamber due to horizontal pressure
    F_SRHorzL = P_HorzLower*TBM.SRShellCSAreaHorz; %Force on lower half of SR Chamber due to horizontal Pressure
    F_SRNorm = F_nTBM*TBM.LSRChamber/TBM.LFront; %Normal Force on the Soil Removal Chamber Shell

    F_fricSRK = (F_SRNorm+F_SRHorzL+F_SRVert+F_SRHorzU)*soilProperties.mu; %Friction force on the soil removal chamber Kinetic
    F_fricSRS = (F_SRNorm+F_SRHorzL+F_SRVert+F_SRHorzU)*soilProperties.mus; %Friction force on the soil removal chamber Static

    %Front Gripper Shell
    F_GShellVert = P_vertUpper *TBM.GripShellCSAreaVert; %Force due to vertical Pressure on gripper shell
    F_GShellHorzU = P_HorzUpper*TBM.GripShellCSAreaHorz; %Force on upper half of Gripper shell due to horizontal pressure
    F_GShellHorzL = P_HorzLower*TBM.GripShellCSAreaHorz; %Force on lower half of Gripper shell due to horizontal Pressure
    F_GShellNorm = F_nTBM*TBM.LGripShell/TBM.LFront; %Normal Force on the gripper Shell

    F_fricGShellK = (F_GShellNorm+F_GShellHorzL+F_GShellHorzU+F_GShellVert)*soilProperties.mu; %Friction force on the Gripper shell Kinetic
    F_fricGShellS = (F_GShellNorm+F_GShellHorzL+F_GShellHorzU+F_GShellVert)*soilProperties.mus; %Friction force on the Gripper shell Static

    %Front Gripper 
    F_GripVert = P_vertUpper *TBM.FrontGripCSAreaVert; %Force due to vertical Pressure on gripper shell
    F_GripHorz = P_HorzMid*TBM.FrontGripCSAreahorz; %Force on upper half of SR Chamber due to horizontal pressure

    F_fricGrip = (F_GripHorz+F_GripVert)*2*TBM.Grippermu; %Friction force on both grippers when retracted

    %Hexapod Shell
    F_HexShellVert = P_vertUpper *TBM.HexapodShellExtCSArea; %Force due to vertical Pressure on hexapod shell
    F_HexShellHorzU = P_HorzUpper*TBM.HexapodShellExtCSArea; %Force on upper half of hexapod shell due to horizontal pressure
    F_HexShellHorzL = P_HorzLower*TBM.HexapodShellExtCSArea; %Force on lower half of hexapod shell due to horizontal Pressure
    F_HexShellNorm = TBM.mHexapod*9.81/2; %Normal Force on the hexapod Shell, only half because other side isn't moving
    %Hexapod Shell
    F_HexShellVertS = P_vertUpper *TBM.HexapodShellContCSArea; %Force due to vertical Pressure on hexapod shell
    F_HexShellHorzUS = P_HorzUpper*TBM.HexapodShellContCSArea; %Force on upper half of hexapod shell due to horizontal pressure
    F_HexShellHorzLS = P_HorzLower*TBM.HexapodShellContCSArea; %Force on lower half of hexapod shell due to horizontal Pressure

    F_fricHexapodK = soilProperties.mu *(F_HexShellVert+F_HexShellHorzU+F_HexShellHorzL+F_HexShellNorm); %Friction Force on Moving Hexapod Part
    F_fricHexapodS = soilProperties.mus *(F_HexShellVertS+F_HexShellHorzUS+F_HexShellHorzLS+F_HexShellNorm);

    Friction.Kinetic = F_fricSRK + F_fricGShellK + F_fricGrip+F_fricHexapodK; %Frictional Force with No Factor of Safety
    Friction.Static = F_fricSRS + F_fricGShellS + F_fricHexapodS; %Max Frictional Force in Static Friction

end