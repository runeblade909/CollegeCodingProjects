%Calculates the friction on the Back of the TBM necessary to be pulled
%forward
function Friction = calcBackFriction(TBM,soilData,soilProperties)
    F_nTBM = TBM.mBack*9.81;%Normal force of the TBM

    % Getting the Necessary Soil Pressures
    P_vertUpper = [zeros(26,1);soilData.sigmaV(1:end-52)]; %Vertical Pressure on the upper half of the TBM Approximating at the centriod of a semicircle
    P_HorzUpper = [zeros(26,1);soilData.sigmaH_o(1:end-52)]; %Horizontal Pressure on upper half of TBM, using at rest pressure because soil is statonary relative to TBM
    P_HorzLower = soilData.sigmaH_o(27:end); %Horizontal Pressure on the bottom of the TBM, using at rest pressure
    P_HorzMid = soilData.sigmaH_o(1:end-26); %Horizontal Pressure at the mid point of the TBM, using at rest pressure

    % BackGripperShell
    F_BackGripShellVert = P_vertUpper *TBM.BackGripCSAreaVert; %Force due to vertical Pressure on back gripper shell
    F_BackGripShellHorzU = P_HorzUpper*TBM.BackGripCSAreaHorz; %Force on upper half of back gripper shell due to horizontal pressure
    F_BackGripShellHorzL = P_HorzLower*TBM.BackGripCSAreaHorz; %Force on lower half of back gripper shell due to horizontal Pressure
    F_BackGripShellNorm = F_nTBM*TBM.LBackGrip/TBM.LBack; %Normal Force on the back gripper shell

    F_fricBackGripS = (F_BackGripShellVert+F_BackGripShellHorzU+F_BackGripShellHorzL+F_BackGripShellNorm)*soilProperties.mus; %Friction force on the back gripper shell Static
    F_fricBackGripK = (F_BackGripShellVert+F_BackGripShellHorzU+F_BackGripShellHorzL+F_BackGripShellNorm)*soilProperties.mu; %Friction force on the back gripper shell kinetic

    % Tunnel Support Storage Section
    F_StorageVert = P_vertUpper *TBM.TSCSArea; %Force due to vertical Pressure on TS Storage
    F_StorageHorzU = P_HorzUpper*TBM.TSCSArea; %Force on upper half of TS Storage due to horizontal pressure
    F_StorageHorzL = P_HorzLower*TBM.TSCSArea; %Force on lower half of TS Storage due to horizontal Pressure
    F_StorageNorm = F_nTBM*(TBM.LStorage+TBM.LReleaseSection)/TBM.LBack; %Normal Force on the TS Storage

    F_fricStorageS = (F_StorageVert+F_StorageHorzU+F_StorageHorzL+F_StorageNorm)*soilProperties.mus; %Friction force on the TS Storage Static
    F_fricStorageK = (F_StorageVert+F_StorageHorzU+F_StorageHorzL+F_StorageNorm)*soilProperties.mu; %Friction force on the TS Storage kinetic

    % Back Gripper
    F_GripVert = P_vertUpper *TBM.BackGripCSAreaVert; %Force due to vertical Pressure on gripper shell
    F_GripHorz = P_HorzMid*TBM.BackGripCSAreaHorz; %Force on upper half of SR Chamber due to horizontal pressure

    F_fricGrip = (F_GripHorz+F_GripVert)*2*TBM.Grippermu; %Friction force on both grippers when retracted

    % Release Arm
    F_ReleaseVert = P_vertUpper *TBM.ReleaseShellCSArea; %Force due to vertical Pressure on release Mechanism
    F_ReleaseHorzU = P_HorzUpper*TBM.ReleaseShellCSArea; %Force on upper half of release Mechanism due to horizontal pressure
    F_ReleaseHorzL = P_HorzLower*TBM.ReleaseShellCSArea; %Force on lower half of release Mechanism due to horizontal Pressure
    
    F_fricReleaseS = (F_ReleaseVert+F_ReleaseHorzU+F_ReleaseHorzL)*TBM.TarpSteelmus; %Friction force on the Release Arm Static
    F_fricReleaseK = (F_StorageVert+F_StorageHorzU+F_StorageHorzL)*TBM.TarpSteelmus-0.01; %Friction force on the Release Arm kinetic

    %Hexapod Shell Static
    F_HexShellVert = P_vertUpper *TBM.HexapodShellExtCSArea; %Force due to vertical Pressure on hexapod shell
    F_HexShellHorzU = P_HorzUpper*TBM.HexapodShellExtCSArea; %Force on upper half of hexapod shell due to horizontal pressure
    F_HexShellHorzL = P_HorzLower*TBM.HexapodShellExtCSArea; %Force on lower half of hexapod shell due to horizontal Pressure
    F_HexShellNorm = TBM.mHexapod*9.81/2; %Normal Force on the hexapod Shell, only half because other side isn't moving
    %Hexapod Shell Kinetic
    F_HexShellVertS = P_vertUpper *TBM.HexapodShellContCSArea; %Force due to vertical Pressure on hexapod shell
    F_HexShellHorzUS = P_HorzUpper*TBM.HexapodShellContCSArea; %Force on upper half of hexapod shell due to horizontal pressure
    F_HexShellHorzLS = P_HorzLower*TBM.HexapodShellContCSArea; %Force on lower half of hexapod shell due to horizontal Pressure
    F_fricHexapodK = soilProperties.mu *(F_HexShellVert+F_HexShellHorzU+F_HexShellHorzL+F_HexShellNorm); %Friction Force on Moving Hexapod Part
    F_fricHexapodS = soilProperties.mus *(F_HexShellVertS+F_HexShellHorzUS+F_HexShellHorzLS+F_HexShellNorm);

    Friction.Kinetic = F_fricBackGripK + F_fricStorageK + F_fricGrip + F_fricReleaseK + F_fricHexapodK; %Kinetic Friction on the entire back section
    Friction.Static = F_fricBackGripS + F_fricStorageS + F_fricReleaseS + F_fricHexapodS; %Static Friction on the entire Back Section


end