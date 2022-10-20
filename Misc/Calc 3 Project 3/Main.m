clc; clear; close all;

syms x y;
r0 = [  1.00,3.50; ...
        4.5,1.5; ...
        2.8,3.7; ...
        2.5,1.8; ...
        3.00,2.00; ...
        1.00,1.00; ...
        4.25,4.10; ... 
        2.75, 3.75 ];
Mountain_Range = mountainRangeInitialize(r0);
%
[xs,ys] = meshgrid(0:.2:5, 0:.2:5);
%Mountain_Eval = double(subs(Mountain_Range, {x,y},{xs,ys}));
Mountain_Gradient = gradient(Mountain_Range);

% Surface Plot
figure(1);
mountainSurfacePlot(Mountain_Range);

% Contour Plot
figure(2);
mountainContourPlot(Mountain_Range,Mountain_Gradient);

% 2D Hike Trail
syms t;
xt = 2.5+1.5*cos((pi/3)*t);
yt = 2.4+1.1*sin((pi/3)*t);
maxT = 6;
Trail_Path = subs(Mountain_Range, {x,y}, {xt,yt});

figure(2)
hold on;
Trail_Line_2D = fplot(xt,yt,[0,maxT],'m');
Trail_Line_2D.LineWidth = 1.5;
hold off;

% 3D Hike Trail
figure(1);
hold on;
Trail_Line_3D = fplot3(xt,yt,Trail_Path,[0,maxT],'m');
Trail_Line_3D.LineWidth = 1.5;
hold off;

% Hour Milestones
Hike_Markers = zeros(floor(maxT),2);
for i = 0:1:maxT
    Hike_Markers(i+1,:) = [subs(xt,t,i),subs(yt,t,i)];
end

% 2D Hour Markers
figure(2);
hold on;
Markers_2D = plot(Hike_Markers(:,1),Hike_Markers(:,2),'.r');
Markers_2D_Start = plot(Hike_Markers(1,1),Hike_Markers(1,2),'.g');
Markers_2D.MarkerSize = 20;
Markers_2D_Start.MarkerSize = 20;
hold off;

% 3D Hour Markers
zMarkers = double(subs(Mountain_Range,{x,y},{Hike_Markers(:,1),Hike_Markers(:,2)}));
figure(1);
hold on;
Markers_3D = plot3(Hike_Markers(:,1),Hike_Markers(:,2),zMarkers,'.r');
Markers_3D_Start = plot3(Hike_Markers(1,1),Hike_Markers(1,2),zMarkers(1),'.g');
Markers_3D.MarkerSize = 20;
Markers_3D_Start.MarkerSize = 20;
hold off;

% 2nd Deriv Test
fxx = diff(diff(Mountain_Range,x),x);
fyy = diff(diff(Mountain_Range,y),y);
fxy = diff(diff(Mountain_Range,x),y);

DSyms = fxx * fyy - (fxy)^2;
D = zeros(length(r0),2);

for i = 1:length(r0)
    D(i,1) = double(subs(DSyms,{x,y},r0(i,:))); % D value
    D(i,2) = double(subs(fxx,{x,y},r0(i,:)));   % fxx value
end

% Peak Values
peaksElevation = zeros(length(r0),1);
for i = 1:length(r0)
    peaksElevation(i) = double(subs(Mountain_Range, {x,y}, r0(i,:)));
end

% Compiled into single variable
criticalPoints(:,[1,2]) = r0;
criticalPoints(:,3) = peaksElevation;
criticalPoints(:,[4,5]) = D;

% Min/Max
maxIndexes = criticalPoints(:,4) > 0 & criticalPoints(:,5) < 0;
minIndexes = criticalPoints(:,4) > 0 & criticalPoints(:,5) > 0;
maxPoints = criticalPoints(maxIndexes,1:3);
minPoints = criticalPoints(minIndexes,1:3);
[maxVal,maxInd] = max(maxPoints(:,3));
[minVal,minInd] = min(minPoints(:,3));
absMax(1,[1,2]) = maxPoints(maxInd,[1,2]);
absMax(1,3) = maxVal;
absMin(1,[1,2]) = minPoints(minInd,[1,2]);
absMin(1,3) = minVal;

% Rates of Change

% Parameter Rates
vt = [diff(xt,t),diff(yt,t),diff(Trail_Path,t)];

% Respect to time
figure(3);
fplot(t,vt(3),[0,maxT]);
title('Change in Elevation w/ Respect to t');
xlim([0,maxT]);
ylim([-4000,4000]);


% Respect to Distance
syms s;
vMagnitude = norm(vt([1,2]));
Distance_Rate = vt(3)/vMagnitude;
figure(4);
fplot(Distance_Rate,[0,maxT]);
xlim([0,6]);
ylim([-2500,2500]);
title('Change in Elevation w/ Respect to distance');

% Max rate of change
time = 0:.01:maxT;
vt_Eval = double(subs(vt(3),t,time));
[absMin,minIndex] = min(vt_Eval);
[absMax,maxIndex] = max(vt_Eval);
changePoints = [absMin,minIndex;absMax,maxIndex];
[absChange,changeIndex] = max(abs(changePoints(:,1)));
maxChangeIni = changePoints(changeIndex,:);
maxChange = changePoints(changeIndex,1);
maxChangeTime = time(maxChangeIni(2));
maxChangePosition = [double(subs(xt,t,maxChangeTime)), ...
                    double(subs(yt,t,maxChangeTime)), ...
                    double(subs(Trail_Path,t,maxChangeTime))];
% Steepest slope
steepestSlopeRate = absMin;
steepestSlopeTime = time(minIndex);
steepestSlopePosition = [double(subs(xt,t,steepestSlopeTime)), ...
                        double(subs(yt,t,steepestSlopeTime)), ...
                        double(subs(Trail_Path,t,steepestSlopeTime))];

                    %{
% Max elevations
hikeEllipse = (x - 2.5)^2/1.5^2+(y - 2.4)^2/1.1^2;

hikeGradient = gradient(hikeEllipse,[x,y]);
%lagrange = cross([Mountain_Gradient;0],[hikeGradient;0]);
lagrange = [Mountain_Gradient(1)*hikeGradient(2),Mountain_Gradient(2)*hikeGradient(1)];
constraintX = solve(hikeEllipse == 1,y);
lagrangeX = solve(lagrange(1) == lagrange(2),y);
%yPoints = solve(lagrangeX == constraintX,y);
%lagrangeX = solve(lagrange == 0,x);
                        %}                        
%
% Time above 10000 ft
elevationIntersection = solve(Trail_Path == 10000,t);

figure(5);
fplot(Trail_Path,[0,maxT]);
hold on;
fplot(10000,[0,maxT]);
hold off;

% Total Distance Walked
%walkDistance = double(int(norm(vt(3)),t,0,maxT));

% Volume of rock inside trail
hikeEllipse = (x - 2.5)^2/1.5^2+(y - 2.4)^2/1.1^2;
hikeEllipseY = solve(hikeEllipse == 1,y);
%rockVolume = double(int(int(Mountain_Range,y,[hikeEllipseY(2),hikeEllipse(1)]),x,[-1.5,1.5]));

% Average height of Rock
%mountainsArea = double(int(int(1,y,[0,5]),x,[0,5]));
%mountainsVolume = double(int(int(Mountain_Range,y,[0,5]),x,[0,5]));
%heightAvg = mountainsVolume/mountainsArea;
