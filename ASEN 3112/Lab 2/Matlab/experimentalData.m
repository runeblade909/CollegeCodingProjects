

clc
close all
clear all

set(0,'defaultTextInterpreter','Latex')
filename = "Distributed_Trial_01";
data = readtable(filename);
data.Properties.VariableNames  = ["LoadingCase","F0","F1","F2","F3D","LVDT"];

%% Displacement
figure;
zeroIdx = data.LoadingCase == 0;
zeroVal = data.LVDT(zeroIdx);
data.displacement = data.LVDT - mean(zeroVal);
scatter(data.LoadingCase,data.displacement)
title("Displacement vs External Load")
xlabel("External Load [lbf]")
ylabel("Displacement [in]")

coeffs = polyfit(data.LoadingCase,data.displacement,3);
x = linspace(0,max(data.LoadingCase),100);
hold on
%plot(x,coeffs(1)*x.^3 + coeffs(2)*x.^2 + coeffs(3)*x + coeffs(4))
%plot(x,coeffs(1)*x.^2 + coeffs(2)*x + coeffs(3))

%Linear Fit/Uncertainty
[slope,interval] = slope_interval(data.LoadingCase,data.displacement);


disp("Displacement:")
disp("    slope = " + slope + " ± " + interval + "(" + -interval*100/slope + "%)");

%% Internal Force
figure;
scatter(data.LoadingCase,data.F3D)
title("Internal Force vs External Load")
xlabel("External Load [lbf]")
ylabel("Internal Force [lbf]")

coeffs = polyfit(data.LoadingCase,data.F3D,2);
x = linspace(0,max(data.LoadingCase),100);
hold on
%plot(x,coeffs(1)*x.^3 + coeffs(2)*x.^2 + coeffs(3)*x + coeffs(4))
%plot(x,coeffs(1)*x.^2 + coeffs(2)*x + coeffs(3))


%Linear Fit/Uncertainty
[slope,interval] = slope_interval(data.LoadingCase,data.F3D);

disp("Internal Force:")
disp("    slope = " + slope + " ± " + interval + "(" + interval*100/slope + "%)");

%% Load Cells
figure;
hold on
scatter(data.LoadingCase,data.F0)
scatter(data.LoadingCase,data.F1)
scatter(data.LoadingCase,data.F2)
scatter(data.LoadingCase,data.F0 + data.F2 + data.F1)
title("Load Cell Forces vs External Load")
xlabel("External Load [lbf]")
ylabel("Force [lbf]")
legend(["F0","F1","F2","\SigmaF"])


%F0
[slope,interval] = slope_interval(data.LoadingCase,data.F0);
disp("F0:")
disp("    slope = " + slope + " ± " + interval + "(" + interval*100/slope + "%)");

%F1
[slope,interval] = slope_interval(data.LoadingCase,data.F1);
disp("F1:")
disp("    slope = " + slope + " ± " + interval + "(" + interval*100/slope + "%)");

%F2
[slope,interval] = slope_interval(data.LoadingCase,data.F2);
disp("F2:")
disp("    slope = " + slope + " ± " + interval + "(" + interval*100/slope + "%)");

%SigmaF
[slope,interval] = slope_interval(data.LoadingCase,data.F2+data.F1+data.F0);
disp("Sigma F:")
disp("    slope = " + slope + " ± " + interval + "(" + interval*100/slope + "%)");


%% Zak Stuff

filename2 = "Ansys Vector data";
data2 = readtable(filename2);

filename3 = "Ansys Vector data (1)";
data3 = readtable(filename3);
data3 = table2array(data3(:,1:3));


% data2.Properties.VariableNames  = ["Node","FX","FY","FZ"];

figure
% scatter(table2array(data2(:,1)),table2array(data2(:,2)))
scatter(table2array(data2(:,1)),table2array(data2(:,3)))
ylim([-.000005,.000005])

figure()

scatter(data3(:,1),data3(:,3)/25.4)
xlabel('Nodes')
ylabel('Displacement[Inch]')
title('Displacement at each Node in Ansys')

%REACTION FORCES ANSYS

W = [-1.1295,55.562,1.0874;0,55.638,1.0040;1.1295,55.638,-1.1580;0,55.562,-0.93343];

figure()
plot(1:4,W(:,2))
xlabel('Supports')
ylabel('Reaction Force [N]')
title('Ansys Support Reactions')


%Bar Forces


filename4 = "Truss3D_Bar_F";
data4 = readtable(filename4);
stuff = table2array(data4(:,1:2))


figure
% scatter(table2array(data2(:,1)),table2array(data2(:,2)))
scatter(stuff(:,1),stuff(:,2))
xlabel('Bar')
ylabel('Bar Force [N]')
title('Ansys Bar Forces')
grid on

function [slope,interval] = slope_interval(x,y)

mdl= fitlm(x,y);
coeff = mdl.Coefficients.Estimate;
uncertainty = coefCI(mdl);
interval = diff(uncertainty(2,:))/2;
slope = coeff(2);
end