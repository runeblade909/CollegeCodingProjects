%driving variables:
maxDiamRing = 0.585; %[m]
minDiamTunnel = 0.5; %[m]
totLengthTunnel = 31.9; %[m]

%known values:
tarpThickness = 0.000254; %[m]
tarpYS = 1.382e7; %[pa]
rodYS = 3.72e8; %[pa]
maxDepthTunnel = 1.8; %[m]

%find:
% rodDiam = 
% distBtwRings = 
% numRings = 

%max sag:
maxSagTop = getMaxSagTop(minDiamTunnel, maxDiamRing);
maxSagBottom = maxSagTop;



function [maxSagTop] = getMaxSagTop(minDiamTunnel,maxDiamRing)
    maxSagTop = (maxDiamRing-minDiamTunnel)/2;
end


