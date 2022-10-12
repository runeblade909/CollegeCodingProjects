clear;
clc;
close all
load('labdata')

R = 287; % j/(kg*k)
A = 1/9.5; % For venturi velocity
T = mean(cell2mat(pitotdata.(namepitot(1))(:,2)));
P = mean(cell2mat(pitotdata.(namepitot(1))(:,1)));
rho = P/(R*T);
%% Pitot VInfinity
% calculating velocity values

% Preallocating 
velVolPitotPT = zeros(60, 2);
velVolVenturiPT = zeros(60, 2);

% Loop through each file
for i = 1:12
    PressureDifferentialPitot = cell2mat(pitotdata(1).(namepitot(i))(:,3)); % getting the Po-P data from pitot data
    PitotVoltage = cell2mat(pitotdata(1).(namepitot(i))(:,7)); % getting the voltage data from pitot data
    
    VenturiPressureDifferential = cell2mat(venturidata(1).(nameventuri(i))(:,3)); % getting the Po-P data from venturi data
    VenturiVoltage = cell2mat(venturidata(1).(nameventuri(i))(:,7)); % getting the voltage data from venturi data

    % Average data by voltage ("collapse into one data point")
    for j = 1:5
        velVolPitotPT(j+(i-1)*5,:) = [mean(PressureDifferentialPitot(500*(j-1)+1:500*j)), mean(PitotVoltage(500*(j-1)+1:500*j))];
        velVolVenturiPT(j+(i-1)*5,:) = [mean(VenturiPressureDifferential(500*(j-1)+1:500*j)), mean(VenturiVoltage(500*(j-1)+1:500*j))];
    end
end



% TODO: convert pressure differentials (col 1) in velVolPitotPT & velVolVentruiPT to velocities
velVolPitotPT(:, 3) = sqrt((2*velVolPitotPT(:, 1)*R*T)/P); % Convert pitot data to velocity
velVolVenturiPT(:, 3) = sqrt((2*velVolVenturiPT(:, 1)*R*T)/(P*(1-A^2))); % Convert venturi data to velocity

%% Find velocity and voltages for manometer data
% man_pitot & man_venturi are columns of velocity and voltage

% Preallocate vectors and iterators for loop
velVolPitotM = NaN(75, 2);
velVolVenturiM = NaN(55, 2);
mpi = 1;
mvi = 1;

% Convert data to cells then pull "names" (venturi or pitot)
manometerdata = table2cell(manometerdata);
mano_names = manometerdata(:, 3);

% Loop through each line in manometerdata
for i = 1:26
    tempname = char(mano_names(i));
        % Determine if current line is venturi or pitot
        if (tempname(1) == 'P')
            % Add all 5 of the height and voltage readings for the line
            velVolPitotM(mpi:mpi+4, 1) = cell2mat(manometerdata(i, 5:2:13)');
            velVolPitotM(mpi:mpi+4, 2) = cell2mat(manometerdata(i, 4:2:12)');
            mpi = mpi + 5;
        else
            % Add all 5 of the height and voltage readings for the line
            velVolVenturiM(mvi:mvi+4, 1) = cell2mat(manometerdata(i, 5:2:13)');
            velVolVenturiM(mvi:mvi+4, 2) = cell2mat(manometerdata(i, 4:2:12)');
            mvi = mvi + 5;
        end
end

velVolPitotM(:, 1) = velVolPitotM(:, 1) * 1000 * 9.81 * 0.0254; % Convert manometer heights to Pa
velVolVenturiM(:, 1) = velVolVenturiM(:, 1) * 1000 * 9.81 * 0.0254; % Convert manometer heights to Pa

velVolPitotM(:, 3) = sqrt((2*velVolPitotM(:, 1)*R*T)/P); % Convert manometer Pa to velocity
velVolVenturiM(:, 3) = sqrt((2*velVolVenturiM(:, 1)*R*T)/(P*(1-A^2))); % Convert manometer Pa to velocity

%% Error Propogation

velVolPitotM = velVolPitotM(velVolPitotM(:,1)~=0, :);
velVolVenturiM = velVolVenturiM(velVolVenturiM(:,1)~=0, :);

PPitotM = polyfit(velVolPitotM(:,2), velVolPitotM(:,3), 1);
PVenturiM = polyfit(velVolVenturiM(:,2), velVolVenturiM(:,3), 1);
PPitotPT = polyfit(velVolPitotPT(:,2), velVolPitotPT(:,3), 1);
PVenturiPT = polyfit(velVolVenturiPT(:,2), velVolVenturiPT(:,3), 1);

% hold on
% plot(PPitotM(1)*velVolPitotM(:,2)+PPitotM(2), velVolPitotM(:,2))
% plot(PPitotPT(1)*velVolPitotPT(:,2)+PPitotPT(2), velVolPitotPT(:,2))
% plot(PVenturiM(1)*velVolVenturiM(:,2)+PVenturiM(2), velVolVenturiM(:,2))
% plot(PVenturiPT(1)*velVolVenturiPT(:,2)+PVenturiPT(2), velVolVenturiPT(:,2))
% legend('PPitotM', 'PPitotPT', 'PVenturiM', 'PVenturiPT')
% hold off

% PPitot = mean([velVolPitotPT(:, 3), velVolPitotM(:, 3)]);
% PVenturi = mean([velVolVenturiPT(:, 3), velVolVenturiM(:, 3)]);

syms Pdiff Tatm Patm
VelocityPitot = sqrt(2.*Pdiff.*((R.*Tatm)./Patm));
VelocityVenturi = sqrt((2*Pdiff(:, 1)*R*Tatm)/(Patm*(1-A)));

dVdPdiffVenturi = diff(VelocityVenturi, Pdiff);
dVdTVenturi = diff(VelocityVenturi, Tatm);
dVdPatmVenturi = diff(VelocityVenturi, Patm);

dVdPdiffPitot = diff(VelocityPitot, Pdiff);
dVdTPitot = diff(VelocityPitot, Tatm);
dVdPatmPitot = diff(VelocityPitot, Patm);


%Given errors from appendix A
dPdiff = .01*6894.745; %[pa]
dTatm = .25; %[C]
dPatm = .015*(250-20)*1000; %[pa]

%Error Ca
dVPitot = sqrt((dVdPdiffPitot*dPdiff)^2+(dVdTPitot*dTatm)^2+(dVdPatmPitot*dPatm)^2);
dVVenturi = sqrt((dVdPdiffVenturi*dPdiff)^2+(dVdTVenturi*dTatm)^2+(dVdPatmVenturi*dPatm)^2);

errorPitotM = subs(dVPitot, {Pdiff, Tatm, Patm}, {velVolPitotM(:,1), T, P});
errorPitotPT = subs(dVPitot, {Pdiff, Tatm, Patm}, {velVolPitotPT(:,1), T, P});
errorVenturiM = subs(dVVenturi, {Pdiff, Tatm, Patm}, {velVolVenturiM(:,1), T, P});
errorVenturiPT = subs(dVVenturi, {Pdiff, Tatm, Patm}, {velVolVenturiPT(:,1), T, P});

figure
hold on
grid on
ylim([-60 70])
errorbar(velVolPitotPT(:,2), PPitotPT(1)*velVolPitotPT(:,2)+PPitotPT(2), errorPitotPT)
errorbar(velVolVenturiPT(:,2), PVenturiPT(1)*velVolVenturiPT(:,2)+PVenturiPT(2), errorVenturiPT)
xlabel('Voltage [V]')
ylabel('Velocity [m/s]')
title('Velocity vs. Voltage - Pressure Transducer')
legend('Pitot', 'Venturi', 'Location','southeast')
hold off

figure
hold on
grid on
ylim([-30 70])
errorbar(velVolPitotM(:,2), PPitotM(1)*velVolPitotM(:,2)+PPitotM(2), errorPitotM)
errorbar(velVolVenturiM(:,2), PVenturiM(1)*velVolVenturiM(:,2)+PVenturiM(2), errorVenturiM)
xlabel('Voltage [V]')
ylabel('Velocity [m/s]')
title('Velocity vs. Voltage - Manometer')
legend('Pitot', 'Venturi', 'Location','southeast')
hold off


%% Process Boundary Layer Data

B = cell(1, 11); % [P diff, ELDy, v] for each port
Bfree = cell(1, 11); % average freestream v for each port

for port = 1:11
    clear Btemp;
    % Loop through each file in each port
    for f = 1:length(nameboundary{port})
        % Save aux P and ELDy data
        Baux = abs(cell2mat(boundarydata{port}.(nameboundary{port}{f})(:, 4)));
        BELDy = cell2mat(boundarydata{port}.(nameboundary{port}{f})(:, 6));
        
        % Average aux & ELDy for every 500 rows
        % Btemp = [aux avg, ELDy avg]
        %Btemp = NaN(length(nameboundary{port})*12, 2);
        for j = 1:12
        Btemp(j+(f-1)*12,:) = [mean(Baux(500*(j-1)+1:500*j)), mean(BELDy(500*(j-1)+1:500*j))];
        end
    end
    % Add avg aux & ELDy data to B for the current port
    B{port} = Btemp;
    % Add 3rd col that is the velocity calculated from the avg aux values (col 1)
    B{port}(:, 3) = sqrt((2.*B{port}(:, 1).*R.*T)./P);
    % average freestream v for each port
    Bfree{port} = mean(B{port}((B{port}(:, 2) > 100),3));
    % Remove rows with freestream data
    B{port}(B{port}(:,2) > 100, :) = []; 
end


%% Plot velocity vs ELD y location for each port

% PB contains the polyfit coefs for each port
% PBvec is PB polyval'ed 
PB = cell(1, 11);
PBvec = cell(1, 11);
figure
hold on
portnames = strings(1,12);
for i = 1:11
    PB{i} = polyfit(B{i}(:, 2), B{i}(:, 3), 3);
    PBvec{i} = polyval(PB{i}, 0:0.1:10);
    plot(0:0.1:10, PBvec{i})
    portnames(i) = strcat('Port ', num2str(i));
end

Bfreeplot = mean(cell2mat(Bfree));
yline(Bfreeplot);
portnames{12} = 'Avg. Freestream';
lgd = legend(portnames, 'Location', 'southeast');
lgd.NumColumns = 3;
xlabel("ELD y location [mm]")
ylabel("Velocity [m/s]")
title("V vs. ELD y")
grid on
hold off


%% Boundary Layer Thickness
% Calculate freestream velocity
Vel_at_B = 0.95 .* cell2mat(Bfree);

% Calculate boundary layer thickness
thickness = zeros(1,11);
for i = 1:11
    tempthickness = roots([PB{i}(1:end-1), PB{i}(end)-Vel_at_B(i)]);
    thickness(i) = tempthickness(end)/1000; %[m]
end
in2m = .0254;
portlocations = [9.05 10.03 11.01 11.99 12.97 13.95 14.93 15.91 16.89 17.87 18.85]; %[in]
portlocations = portlocations*in2m;

% Best fit
PBThickness = polyfit(portlocations, thickness, 2);
PBThicknessvec = polyval(PBThickness, portlocations);

% Plot Thickness vs Port locations
figure
hold on
x_lower = portlocations(1)*1000;
x_higher = portlocations(end)*1000;
xlim([x_lower x_higher])
a = area(portlocations*1000, PBThicknessvec*1000, 'FaceColor','b', 'EdgeColor','none');
a.FaceAlpha = 0.3;
scatter(portlocations*1000, thickness*1000, 'r*');
xlabel("Port Location [mm]")
ylabel("Boundary Layer Thickness [mm]")
title("Boundary Layer vs. Streamwise Coordinate")
hold off


%% Turbulent and Laminar Boundary Thicknesses
% Compares the experimentally determined boundary layer thickness
% with the calculated turbulent and laminar boundary thicknesses
Vinfinity = cell2mat(Bfree);
mu = 1.7894 * 10^(-5); % [kg/m*s] viscosity at sea level
Rex = (rho .* Vinfinity .* portlocations)/mu; % Local reynolds number

% Calculate turbulent and laminar thicknesses
turbulent = Rex;
laminar = Rex;
turbulent = (0.37 .* portlocations) ./ turbulent.^0.2;
laminar = (5.2 .* portlocations) ./ laminar.^0.5;

ThicknessP = polyfit(portlocations, thickness, 1);
Thicknessvec = polyval(ThicknessP, portlocations);
% Plot thickness vs port location
figure
hold on
xlim([x_lower x_higher])
plot(portlocations*1000, Thicknessvec*1000, 'r', 'LineWidth',1.5);
plot(portlocations*1000, turbulent*1000, 'g', 'LineWidth',1.5);
plot(portlocations*1000, laminar*1000, 'b', 'LineWidth',1.5);
scatter(portlocations*1000, thickness*1000, 'r*');
xlabel("Port Location [mm]")
ylabel("Boundary Layer Thickness [mm]")
title("Comparing Possible Boundary Layers")
legend('Results', 'Turbulent', 'Laminar', 'location', 'northwest')
hold off

% Based on the graph the flow in the boundary layer is laminar

%close all; warning('close all at the beginning of part 7!');


%% Part 7
%% Calculating trailing edge pressure port value
% Get rid of 302_3 because it has one less row
airfoildata = rmfield(airfoildata, nameairfoil(11));
nameairfoil(11) = [];

n = length(nameairfoil);
AP_PTD = zeros(n,3);    % AirfoilPressure_Pressure,Temp,Density
AP = zeros(n, 20);
% Loop through each file
for i = 1:n
    % All scanivalve & AOA data per file
    tempdata = cell2mat(airfoildata.(nameairfoil(i))(:, 4:23));
    tempPTD = cell2mat(airfoildata.(nameairfoil(i))(:, 1:3));
    AP_PTD(i,:) = mean(tempPTD);
    % Average data by AOA ("collapse into one data point")
    for j = 1:12
        AP(j+(i-1)*12,:) = mean(tempdata(20*(j-1)+1:20*j, :));
    end
end
AP = sortrows(AP, 1);
AP_PTD = mean(AP_PTD);

% Separate AOA based on airspeed
AP9 = AP(abs(AP(:,1)-9)<=2, :);
AP16 = AP(abs(AP(:,1)-16)<=2, :);
AP34 = AP(abs(AP(:,1)-34)<=2, :);

% Sort AOA matrices by angle of attack
AP9 = sortrows(AP9,20);
AP16 = sortrows(AP16,20);
AP34 = sortrows(AP34,20);

% Linearly fit ports 8,10 & 12,14 / scanivalve ports #s 8,9 & 10,11 / columns 11-14
% Using the linear fit, we can find out trailing edge pressure value

top16p = zeros(100,2);
top34p = zeros(100,2);
bottom16p = zeros(100,2);
bottom34p = zeros(100,2);
for i = 1:100
    % Last two ports on the top
    top16p(i,:) = polyfit([2.1 2.8],AP16(i,[11 12]),1);
    top34p(i,:) = polyfit([2.1 2.8],AP34(i,[11 12]),1);
    
    % Last two ports on the bottom
    bottom16p(i,:) = polyfit([2.1 2.8],AP16(i,[14 13]),1);
    bottom34p(i,:) = polyfit([2.1 2.8],AP34(i,[14 13]),1);
end

top9p = zeros(100,2);
bottom9p = zeros(100,2);
for i = 1:100
    % Last two ports on the top
    top9p(i,:) = polyfit([2.1 2.8],AP9(i,[11 12]),1);
    % Last two ports on the bottom
    bottom9p(i,:) = polyfit([2.1 2.8],AP9(i,[14 13]),1);   
end

% Calculate the pressure at the trailing edge (3.5 inch chord) based on the top and bottom
% linear fits
trailPtop9 = top9p(:,1).*3.5+top9p(:,2);
trailPtop16 = top16p(:,1).*3.5+top16p(:,2);
trailPtop34 = top34p(:,1).*3.5+top34p(:,2);

trailPbottom9 = bottom9p(:,1).*3.5+bottom9p(:,2);
trailPbottom16 = bottom16p(:,1).*3.5+bottom16p(:,2);
trailPbottom34 = bottom34p(:,1).*3.5+bottom34p(:,2);

% Average the trailing edge pressures from the top and bottom fits
TrailingEdgeP9 = mean([trailPtop9 trailPbottom9], 2);
TrailingEdgeP16 = mean([trailPtop16 trailPbottom16], 2);
TrailingEdgeP34 = mean([trailPtop34 trailPbottom34], 2);

% Add trailing edge pressure value AP
AP9 =  [AP9(:, 1:12) TrailingEdgeP9 AP9(:, 13:20)];
AP16 = [AP16(:, 1:12) TrailingEdgeP16 AP16(:,13:20)];
AP34 = [AP34(:, 1:12) TrailingEdgeP34 AP34(:,13:20)];
% AP looks like:
% airspeed, pitot dynamic, aux dynamic, scanivalves, trailing edge,  AOA

% Collapse AP data by AOA
% https://www.mathworks.com/matlabcentral/answers/354507-how-do-i-average-data-in-a-matrix-that-has-the-same-coordinates-as-a-different-data-point
[u9,~,j9] = unique(AP9(:,21),'rows','stable');
[u16,~,j16] = unique(AP16(:,21),'rows','stable');
[u34,~,j34] = unique(AP34(:,21),'rows','stable');

AP9new = zeros(length(u9), 21);
AP16new = zeros(length(u16), 21);
AP34new = zeros(length(u34), 21);

AP9new(:,21) = u9;
AP16new(:,21) = u16;
AP34new(:,21) = u34;

% Average each column of rows that have the same AOA
for i = 1:20
    AP9new(:,i) = accumarray(j9,AP9(:,i),[],@mean);
    AP16new(:,i) = accumarray(j16,AP16(:,i),[],@mean);
    AP34new(:,i) = accumarray(j34,AP34(:,i),[],@mean);
end

AP9 = AP9new;
AP16 = AP16new;
AP34 = AP34new;
% Clean up workspace
clear AP9new AP16new AP34new result9 result 16 result 34 u9 u16 u34 j9 j16 j34;


%% Calculate Coef. of Pressure at each port
Cp9 = zeros(32,17);
Cp16 = zeros(32,17);
Cp34 = zeros(32,17);

for i = 1:17
    Cp9(:,i) = (AP9(:,i+3))./(0.5 .* AP_PTD(3) .* AP9(:,1).^2);
    Cp16(:,i) = (AP16(:,i+3))./(0.5 .* AP_PTD(3) .* AP16(:,1).^2);
    Cp34(:,i) = (AP34(:,i+3))./(0.5 .* AP_PTD(3) .* AP34(:,1).^2);
end

% Portlocations = [x;y] for scanivalve 1-16 and trailing edge
portlocations = [0 0.175 0.35 0.7 1.05 1.4 1.75 2.1 2.8 3.5 2.8 2.1 1.4 1.05 0.7 0.35 0.175 ;
                 0.14665 0.33075 0.4018 0.476 0.49 0.4774 0.4403 0.38325 0.21875 0 0 0 0 0 0.0014 0.0175 0.03885];
%percent chord
xc = portlocations(1,:) ./ 3.5;

% close all

figure %9[m/s]
grid on
hold on

plot(xc([1:17 1]), Cp9(6,[1:17 1]),'m-')
plot(xc([1:17 1]), Cp9(11,[1:17 1]),'c-')
plot(xc([1:17 1]), Cp9(16,[1:17 1]),'r-')
plot(xc([1:17 1]), Cp9(21,[1:17 1]),'b-')
plot(xc([1:17 1]), Cp9(26,[1:17 1]),'g-')

axis ij

legend('-10','-5','0','5','10','Location','southeast')
hold off
title('Coeff. of Pressure (9[m/s]) vs. % Chord')
xlabel('% of Chord length')
ylabel('Coeff. of Pressure')

figure %16[m/s]
grid on
hold on

plot(xc([1:17 1]), Cp16(6,[1:17 1]),'m-')
plot(xc([1:17 1]), Cp16(11,[1:17 1]),'c-')
plot(xc([1:17 1]), Cp16(16,[1:17 1]),'r-')
plot(xc([1:17 1]), Cp16(21,[1:17 1]),'b-')
plot(xc([1:17 1]), Cp16(26,[1:17 1]),'g-')

axis ij

legend('-10','-5','0','5','10','Location','northeast')
hold off
title('Coeff. of Pressure (16 [m/s]) vs. % Chord')
xlabel('% of Chord length')
ylabel('Coeff. of Pressure')

figure %34[m/s]
grid on
hold on

plot(xc([1:17 1]), Cp34(6,[1:17 1]),'m-')
plot(xc([1:17 1]), Cp34(11,[1:17 1]),'c-')
plot(xc([1:17 1]), Cp34(16,[1:17 1]),'r-')
plot(xc([1:17 1]), Cp34(21,[1:17 1]),'b-')
plot(xc([1:17 1]), Cp34(26,[1:17 1]),'g-')

axis ij

legend('-10','-5','0','5','10','Location','northeast')
hold off
title('Coeff. of Pressure (34[m/s]) vs. % Chord')
xlabel('% of Chord length')
ylabel('Coeff. of Pressure')

figure
grid on
hold on

plot(xc([1:17 1]), Cp9(21,[1:17 1]),'r-')
plot(xc([1:17 1]), Cp16(21,[1:17 1]),'b-')
plot(xc([1:17 1]), Cp34(21,[1:17 1]),'g-')

axis ij

legend('9 [m/s]','16[m/s]','34[m/s]','Location','northeast')
hold off
title('(5 degree AOA) Coeff. of Pressure vs. % Chord')
xlabel('% of Chord length')
ylabel('Coeff. of Pressure')



%% Coeff. of N and A

Cn9 = zeros(1,32);
Ca9 = zeros(1,32);
Cn16 = zeros(1,32);
Ca16 = zeros(1,32);
Cn34 = zeros(1,32);
Ca34 = zeros(1,32);

%Calculate Coeff. Normal and Coeff. Axial (Pg.20 Eq.16-17)
for j =1:32
    for i = 1:16
    DeltaX = portlocations(1,i+1) - portlocations(1,i);
    DeltaY = portlocations(2,i+1) - portlocations(2,i);
    
    Cn9(j) = Cn9(j) + (1/2)*(Cp9(j,i) + Cp9(j,i+1)) *(DeltaX/3.5);
    Ca9(j) = Ca9(j) + (1/2)*(Cp9(j,i) + Cp9(j,i+1)) *(DeltaY/3.5);

    Cn16(j) = Cn16(j) + (1/2)*(Cp16(j,i) + Cp16(j,i+1)) *(DeltaX/3.5);
    Cn34(j) = Cn34(j) + (1/2)*(Cp34(j,i) + Cp34(j,i+1)) *(DeltaX/3.5);

    Ca16(j) = Ca16(j) + (1/2)*(Cp16(j,i) + Cp16(j,i+1)) *(DeltaY/3.5);
    Ca34(j) = Ca34(j) + (1/2)*(Cp34(j,i) + Cp34(j,i+1)) *(DeltaY/3.5);
    
    end
end
%Make sum of Cn Negative like it is in the formula
Cn9 = -Cn9;
Cn16 = -Cn16;
Cn34 = -Cn34;

Cd9 = zeros(1,32);
Cl9 = zeros(1,32);
Cd16 = zeros(1,32);
Cl16 = zeros(1,32);
Cd34 = zeros(1,32);
Cl34 = zeros(1,32);

%% Coeff. of Lift and Coeff. of Drag Calculations
alpha = AP9(:,21);
for i = 1:32
    Cl9(i) = Cn9(i).*cosd(alpha(i)) - Ca9(i).*sind(alpha(i));
    Cd9(i) = Cn9(i).*sind(alpha(i)) + Ca9(i).*cosd(alpha(i));

    Cl16(i) = Cn16(i).*cosd(alpha(i)) - Ca16(i).*sind(alpha(i));
    Cl34(i) = Cn34(i).*cosd(alpha(i)) - Ca34(i).*sind(alpha(i));
    Cd16(i) = Cn16(i).*sind(alpha(i)) + Ca16(i).*cosd(alpha(i));
    Cd34(i) = Cn34(i).*sind(alpha(i)) + Ca34(i).*cosd(alpha(i));
end
% Plot Cl vs AOA & Cd vs AOA

figure
hold on
grid on
plot(alpha, Cl9);
plot(alpha, Cl16);
plot(alpha, Cl34);
legend('9 m/s', '16 m/s', '34 m/s', 'location', 'northwest');
title('Coef. of Lift vs. Angle of Attack');
xlabel('Angle of Attack [deg]');
ylabel('Coefficent of Lift');
xlim([-15 16])


figure
hold on
grid on
plot(alpha, Cd9);
plot(alpha, Cd16);
plot(alpha, Cd34);
legend('9 m/s', '16 m/s', '34 m/s', 'location', 'northwest');
title('Coef. of Drag vs. Angle of Attack');
xlabel('Angle of Attack [deg]');
ylabel('Coefficent of Drag');
xlim([-15 16])



