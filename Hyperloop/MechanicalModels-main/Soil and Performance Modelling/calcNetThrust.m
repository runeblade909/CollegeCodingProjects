% Calculates the thrust needed for different soil conditions,penetration
% depths ranging from 1 for full to 0 for no pentetration, and the geometric 
% properties of the cutting head, Computes across all of the depths for
% which the soil data has been computed
function Thrust = calcNetThrust(cuttingHead,soilData,soilProperties,Penetration)
% Inputs: 
% CuttingHead: A struct containing the geometry of the cutting head
% soilData: A struct containing the computed stengths of the soil in various modes and planes
% Soil Properties: A Struct containing properties of the soil that do not
% vary
% Penetration, The Penetration of the Tooth, can range from 0-1

% Outputs:
% Thrust: The Net thrust required to penetrate the cutting head

        F_TeethInner = soilData.sigmaInnerTooth.*cuttingHead.TeethArea*Penetration + cuttingHead.TeethSA*Penetration.^2*soilData.sigmaH_p*soilProperties.mu; %force on each inner tooth
        F_TeethOuter = soilData.sigmaOuterTooth.*cuttingHead.TeethArea*Penetration + cuttingHead.OuterTeethSA*Penetration.^2*soilData.sigmaH_p*soilProperties.mu; %force on each outer tooth
        F_Spike = soilData.sigmaSpike.*(cuttingHead.SpikeRad-(cuttingHead.TeethPenetration*(1-Penetration))).^2*pi + cos(cuttingHead.SpikeAngle)*cuttingHead.SpikeSA*soilData.sigmaH_p*soilProperties.mu*Penetration.^2; %Penetration force required to push the spike into the dirt
        F_HeadPlate = soilData.sigmaH_p.*cuttingHead.A_HeadBare*Penetration.^2;

    F_Teeth = F_TeethInner * cuttingHead.InnerTeethNum + F_TeethOuter*cuttingHead.OuterTeethNum + F_Spike;
    Thrust = (F_HeadPlate + F_Teeth);
    %Thrust.Design = Thrust.Total * FS;
    
 end