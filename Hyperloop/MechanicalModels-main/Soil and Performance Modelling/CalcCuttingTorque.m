%Calculates the necessary torque given soil conditions and forward thrust
function Torque = CalcCuttingTorque(cuttingHead,SoilData,SoilProperties,Penetration,FS)

    if Penetration == 1 % Full Penetration
        F_fHead = SoilData.mu*SoilProperties.sigmaH_p*cuttingHead.A_HeadBare; %Force of friction on the cutting head plate
        F_Tooth_Rotation = SoilProperties.sigmaInnerTooth * cuttingHead.InnerTeethCSArea; %Force required to pierce the soil for a single inner tooth
        F_OuterTooth_Rotation = SoilProperties.sigmaOuterTooth * cuttingHead.OuterTeethCSArea; %Force required to pierce the soil for a single outer tooth
        F_fInnerTeeth = cuttingHead.TeethSA * SoilProperties.sigmaH_p * SoilData.mu; %Force of friction on a single inner tooth
        F_fOuterTeeth = cuttingHead.OuterTeethSA * SoilProperties.sigmaH_p * SoilData.mu; %Force of friction on a single inner tooth
        F_fSpike = sin(cuttingHead.SpikeAngle)*cuttingHead.SpikeSA*SoilProperties.sigmaH_p*SoilData.mu; %Force of friction on the spike
    else % .5<Pen<1
        if 1 > Penetration && Penetration>= 0.5
            CSArea = 1/3 * cuttingHead.InnerTeethCSArea + (Penetration-0.5) * 2/3 * cuttingHead.InnerTeethCSArea;
        else %Pen<0.5
            CSArea = 1/3 * cuttingHead.InnerTeethCSArea * (Penetration * 2)^2;
        end
        
        %Changed for pen variation
        
        OuterCSArea = Penetration * cuttingHead.OuterTeethCSArea; %Effective Cross Sectional area of outer tooth
        F_fHead = SoilData.mu*SoilProperties.sigmaH_p*cuttingHead.A_HeadBare* Penetration^2;
        F_Tooth_Rotation = SoilProperties.sigmaInnerTooth * CSArea;
        F_OuterTooth_Rotation = SoilProperties.sigmaOuterTooth * OuterCSArea;
        F_fInnerTeeth = cuttingHead.TeethSA *Penetration^2* SoilProperties.sigmaH_p * SoilData.mu;
        F_fOuterTeeth = cuttingHead.OuterTeethSA* Penetration^2 * SoilProperties.sigmaH_p * SoilData.mu;
        F_fSpike = sin(cuttingHead.SpikeAngle)*cuttingHead.SpikeSA*SoilProperties.sigmaH_p*SoilData.mu * Penetration^2;
    end
    
    %TorqueOuterEdge = SoilData.mu*SoilProperties.sigmaV*cuttingHead.d*cuttingHead.Thickness;
    F_ToothInnerTotal = F_fInnerTeeth+F_Tooth_Rotation; %Total force required to rotate a single inner tooth
    F_OuterToothTotal = F_fOuterTeeth+F_OuterTooth_Rotation; %Total Force required to rotate a single outer tooth
    TorqueHead = F_fHead*cuttingHead.Midpoint; %Torque due to friction on the cutting head
    TeethTorqueMultiplier = sum(cuttingHead.TeethLocationNum.*cuttingHead.TeethDist); % AMount by which to multiply the force on a single inner toothe to get the torque on all of the inner teeth
    TorqueInnerTeeth = F_ToothInnerTotal*TeethTorqueMultiplier; %Total torque required for all of the inner teeth
    TorqueOuterTeeth = F_OuterToothTotal*cuttingHead.d/2* cuttingHead.OuterTeethNum; %Total torque required for all of the outer teeth
    TorqueSpike = F_fSpike * sqrt(cuttingHead.SpikeRad^2/2); %Torque due to the friction on the spike
    TorqueTeeth = TorqueInnerTeeth + TorqueOuterTeeth + TorqueSpike; %Total torque due to all of the teeth and the spike
    % TorqueTotal = (TorqueTeeth + Torque.Head);
    % Torque.NonRecessed = TorqueTotal + TorqueOuterEdge;
    Torque = (TorqueTeeth + TorqueHead) * FS; %Total torque required with the desired Factor of Safety
end