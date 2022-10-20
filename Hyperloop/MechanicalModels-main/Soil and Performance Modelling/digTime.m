% Function to calculate the time to completion for dig given the forces and
% thrusts of the dig

function [Time,VFront] = digTime(ForceFront,ForceBack,path,Actuator)
    ForceFrontMidPanel = ForceFront(1:end-1) + diff(ForceFront)/2; % Finding the thrust force at the mid point between each path panel
    ForceBackMidPanel = ForceBack(1:end-1) + diff(ForceBack)/2; % Finding the thrust force at the mid point between each path panel
    VFront = Actuator.maxSpeed + ForceFrontMidPanel*Actuator.SpeedConstant; %Speed during Digging Calculated via the Hexapod
    VBack = Actuator.maxSpeed + ForceBackMidPanel*Actuator.SpeedConstant; %Speed during pulling Calculated via the Hexapod
    diffDist = sqrt(sum(diff(path).^2,2));
    Time.ActiveDigTime = sum(diffDist./VFront); %Active Dig Time when the cutting head is actively diggint min
    Time.PullTime = sum(diffDist./VBack); %Pull time when the back is being pulled forward min
    numStrokes = sum(diffDist)/Actuator.StrokeLength;
    Time.Sample = numStrokes*Actuator.sampleTime/60 - Time.PullTime; %Time required to be stainary for the barometers
    Time.GripperExtendTime = Actuator.gripperStokeTime/60 * numStrokes; %Time the grippers spent extending
    Time.Stationary = (1/Actuator.DutyCycle -1)*(Time.PullTime+Time.ActiveDigTime)-Time.GripperExtendTime;
    if Time.Stationary <0
        Time.Stationary = 0;
    end
    Time.Sample = numStrokes*Actuator.sampleTime/60 - Time.PullTime - Time.Stationary; %Time required to be stainary for the barometers
    if Time.Sample < 0
        Time.Sample = 0;
    end

    Time.Stationary = Time.Stationary + Time.Sample;

    Time.TotalTime = Time.ActiveDigTime + Time.PullTime + Time.GripperExtendTime+Time.Stationary; %Total time to completion of the Dig

end