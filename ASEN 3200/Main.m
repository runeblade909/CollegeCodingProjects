%Main
%Christopher Lolkema
%ASEN3200 Orbital Projects Part 1 Code
tic %Start Time Elapsed Calc
%% Set Up Coding Environment
clear;
clc;
close all;

%% Constants
mu = 3.98e14/1000^3; %Gravitational Parameter [km^3/sec^2]
J2 = 1.087e-3; %Pertibations [dimensionless]
Re = 6378; %Radius of Earth [km]
timeVec = [0:60:86400];
timeLength = length(timeVec);

%% Load World Cities Data and CoastLine Data
cities = readtable("worldcities.csv");
coastline = readtable("world_coastline_low.txt");

%% Limit Amount of Cities for Testing
cities = cities(1:20:41000,:);

%% Load JSON File
filename = "Constellation.json";
[num_launches, num_spacecraft, satellite_list] = loadConstellation(filename);

fprintf('Constellation Selected.\n\n')

%% Propagate State
hold on
for j = 1:num_spacecraft
    t_0 = 0;
    k = 1;
    orbits = zeros(timeLength,6);
    for i = timeVec
        t = i;
        orbits(k,:) = propagateState(satellite_list(j).oe0,t,t_0,mu,J2,Re);
        k = k + 1;
    end
    satellite_list(j).orbitmat = orbits;
end

fprintf('Trajectories Calculated.\n\n')

%% Convert ECEF
for i = 1:num_spacecraft
    statevec = satellite_list(i).orbitmat;
    [azimuth,elevation,r] = cart2sph(statevec(:,1),statevec(:,2),statevec(:,3));
    azimuth = (180/pi * azimuth) - (360/86400).*t + 180/pi * satellite_list(i).AzimuthShift;
    elevation = 180/pi * elevation;
    [maxr,maxi] = max(azimuth);
    [minr,mini] = min(azimuth);
    while maxr > 180 || minr < -180
        for j = 1:length(azimuth)
            if azimuth(j) > 180
                azimuth(j) = azimuth(j) - 360;
            end
            if azimuth(j) < -180
                azimuth(j) = azimuth(j) + 360;
            end
        end
        [maxr,maxi] = max(azimuth);
        [minr,mini] = min(azimuth);
    end
    [x,y,z] = sph2cart(pi/180*azimuth, pi/180*elevation, r);
    satellite_list(i).ecef  = [x,y,z];
end

fprintf('Coordinate Frame Converted\n\n')

cost = num_launches * 100 + num_spacecraft * 5;


%% Compute LoS
figure(1)
hold on
[numCities,~] = size(cities);
%Plot Cities
[siteRX,siteRY,siteRZ] = sph2cart(pi/180*cities{:,4},pi/180*cities{:,3},Re);
scatter3(siteRX,siteRY,siteRZ,.25,"*r")

angleTest = 15;
angleTest = sind(15);
numOrbitalPoints = timeLength;
cityStruct = struct([]);




for i = 1:numCities
    siteR = [siteRX(i),siteRY(i),siteRZ(i)];
    sumVec = zeros(numOrbitalPoints,1);
    for j = 1:numOrbitalPoints
        sumVal = 0;
        for k = 1:num_spacecraft
            sumVal = sumVal + testLoS(siteR,satellite_list(k).ecef(j,:),angleTest);
        end
        sumVec(j) = sumVal;
    end
    cityStruct(i).LoS = sumVec;
end

fprintf('Line of Sight Calculated for all cities.\n\n')

%% Plot Data
view(3)
%Plot Earth Sphere
[X,Y,Z] = sphere();
X = Re*X;
Y = Re*Y;
Z = Re*Z;
mesh(X,Y,Z)
[earthX,earthY,earthZ] = sph2cart(pi/180 * coastline{:,1},pi/180 * coastline{:,2},Re);
scatter3(earthX,earthY,earthZ,1,'.k')
%Plot Orbits
for i = 1:num_spacecraft
    %plot3(satellite_list(i).ecef(:,1),satellite_list(i).ecef(:,2),satellite_list(i).ecef(:,3))
    %plot3(satellite_list(i).orbitmat(:,1),satellite_list(i).orbitmat(:,2),satellite_list(i).orbitmat(:,3))
    %scatter3(satellite_list(i).ecef(1,1),satellite_list(i).ecef(1,2),satellite_list(i).ecef(1,3))
end

hold off
%Set square limits
axis equal

%% Plot City LoS Graph
% figure(2)
% hold on
% for i = 1:numCities
%     plot(timeVec,cityStruct(i).LoS)
% end
% ylim([0,num_spacecraft+1])
% hold off
% xlabel("Time [sec]")
% ylabel("Satellites with Line of Sight")
% title("Plotting Line of Sight with Cities as a Function of Time")

%% Ground Track Plot

%Plot Earth
% figure(2)
% plot(coastline{:,1},coastline{:,2})
% xlim([-180,180])
% ylim([-90,90])
%Plot Orbits
colorVec = linspace(1,100,1);
colorVec = zeros(1,numCities);
figure
for j = 1:length(satellite_list(i).ecef(:,1))
    plot(coastline{:,1},coastline{:,2})
    xlim([-180,180])
    ylim([-90,90])
    hold on
    for k = 1:numCities
        colorVec(k) = cityStruct(k).LoS(j);
    end
    scatter(cities{:,4},cities{:,3},5,colorVec,'filled')
    for i = 1:num_spacecraft
        [Az,Ele,~] = cart2sph(satellite_list(i).ecef(j,1),satellite_list(i).ecef(j,2),satellite_list(i).ecef(j,3));
        scatter(180/pi*Az,180/pi*Ele,'.r')
    end
    hold off
    print(['Frame ' num2str(j)], '-dpng', '-r150');
   
end

%% Export GIF
GifName = 'SatelliteConstellation.gif';
delay = 0.01;    % Delay between frames (s)
for ii = 1:length(satellite_list(i).ecef(:,1))
    [A, ~] = imread(['Frame ' num2str(ii) '.png']);
    [X, map] = rgb2ind(A, 256);
    if ii == 1
        imwrite(X, map, GifName, 'gif', 'LoopCount', inf, 'DelayTime', delay)
    else
        imwrite(X, map, GifName, 'gif', 'WriteMode', 'append', 'DelayTime', delay)
    end
end


%% Time Elapsed
timeElapse = toc/60;
disp(num2str(timeElapse))
disp(' Minutes Elapsed')

