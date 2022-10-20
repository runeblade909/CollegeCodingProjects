%Calculating the Max Teeth Penetration at any given Depth
function Penetration = calcPenetration(cuttingHead,Thrust,Actuator,Penetrations,depth,path)
% Inputs:
% cuttingHead: The Geomatric properties of the cutting head
% Thrust: The Thrust required for varying Penetrations and depths
% Actuator: Actuator Data with speed and force calibrations
% Penetrations: Penetrations of the teeth for a given force value
% depth: Depth for a given thrust value
% path that the TBM is following

    path(:,3) = abs(path(:,3)); %Making sure all of the depths are positive
    % Solving for a required Penetration at each point in the path
    Penetration = zeros(length(path(:,3)),1);
    for i = 1:length(path(:,3))
        DepthIdx = find(depth > path(i,3),1,"first"); %Index at which to check for the thrust
        ThrustDepth = Thrust(DepthIdx,:); %Thrust for a given depth and all penetrations
        Penetration(i) = findConvergence(ThrustDepth,cuttingHead.RPM,cuttingHead.TeethPenetration,Penetrations,Actuator); %Finding the max penetration for a given depth
    end

    %Finds the convergence of the penetration for a given depth
    function PenetrationMax = findConvergence(Thrust,RPM,ToothLength,pen,Actuator)
        % Inputs:
        % Thrust: Thrust for a given depth for all penetrations
        % RPM: RPM of the Cutting Head
        % ToothLength: Cutting length of a tooth on the cutting head
        % pen: Possible Penetrations
        % Actuator: Constants from the actuator datasheet

        % Outputs:
        % Max penetration for a given depth

        V_cuttingHead = pen*ToothLength*RPM; %Speed calculated via the Cutting Head m/min
        V_Hexapod = Actuator.maxSpeed + Thrust'*Actuator.SpeedConstant; %Speed Calculated via the Hexapod
        V_net = V_Hexapod-V_cuttingHead; %Difference in calculation between the two methods
        PenetrationMax = pen(find(V_net>0,1,"last")); %Finding the highest penetration where the hexapod can keep up with the cutting head
    end
end