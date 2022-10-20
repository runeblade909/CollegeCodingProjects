% Function to Find the soil properties at each depth for a given soil
function [ResultSoilData,coefficients] = ModelSoil(soilData,cuttingHead)
    %Terzaghi Soil Model
    ResultSoilData.depths = (0:0.005:2.5)'; %Depths for which to compute the soil properties
    coefficients.Ko = 1-sind(soilData.phi); %At rest lateral pressure coefficient
    coefficients.Ka = (1-sind(soilData.phi))/(1+sind(soilData.phi)); %Active pressure coefficient, soil moving away
    coefficients.Kp = (1+sind(soilData.phi))/(1-sind(soilData.phi)); %Passive Pressure coefficient, soil in compression
    B = cuttingHead.d .* (1 + tand(45 - soilData.phi/2)); % Silo width, using GB model for this [m]
    ResultSoilData.sigmaV = B * soilData.gamma * ((1/(coefficients.Ka*tand(soilData.phi))) .* (1-exp(-coefficients.Ka.*(ResultSoilData.depths./B)*tand(soilData.phi)))); %Vertical Soil Pressure from Terzaghis soil model, Equation 6b on Page 72 of terzaghi book
    ResultSoilData.sigmaH_p = coefficients.Kp*ResultSoilData.sigmaV + 2*soilData.c*sqrt(coefficients.Kp); %Calculating the passive lateral pressure using Rankine Earth pressure coefficient
    ResultSoilData.sigmaH_a = coefficients.Ka*ResultSoilData.sigmaV - 2*soilData.c*sqrt(coefficients.Ka); %Calculating the active lateral earth pressure using Rankine's earth pressure coefficient
    ResultSoilData.sigmaH_o = coefficients.Ko*ResultSoilData.sigmaV; %Calculating the at rest soil pressure using Rankines Earth Pressure coefficients

    %Strength on the planes of various parts of the cutting head
    ResultSoilData.sigmaPlate = soilData.c+ResultSoilData.sigmaH_p*tand(90-soilData.phi); %normal to the main cutting head plate
    ResultSoilData.sigmaInnerTooth = soilData.c+(sqrt((sind(cuttingHead.TeethAngle)*ResultSoilData.sigmaH_p).^2+(cosd(cuttingHead.TeethAngle)*ResultSoilData.sigmaV).^2))*tand(90-soilData.phi); %Soil Strength Normal to the blade of the inner teeth
    ResultSoilData.sigmaOuterTooth = soilData.c+(sqrt((sind(cuttingHead.OuterTeethAngle)*ResultSoilData.sigmaH_p).^2+(cosd(cuttingHead.OuterTeethAngle)*ResultSoilData.sigmaV).^2))*tand(90-soilData.phi); %Soil Strength Normal to the blade of the outer teeth
    ResultSoilData.sigmaSpike = soilData.c+(sqrt((sind(cuttingHead.SpikeAngle)*ResultSoilData.sigmaH_p).^2+(cosd(cuttingHead.SpikeAngle)*ResultSoilData.sigmaV).^2))*tand(90-soilData.phi); %Soil Strength Normal to the spike
    
end