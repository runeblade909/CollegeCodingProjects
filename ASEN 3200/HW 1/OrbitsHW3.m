clc
clear all
close all

Earth.r = 6378; %km
Earth.m = 5.97219*10^24;%kg
Earth.G =  6.67259*10^-11;
mu = Earth.m*Earth.G;


x = load('catalog_TLEs.mat');%Gives Cell Array of Obj

y = x.catalog_TLEs{1};

for i = 1:length(x.catalog_TLEs)
Eccentricity(i) = x.catalog_TLEs{i}.eccentricity;       %Grabs Struct of Obj
Semimajor(i) = x.catalog_TLEs{i}.semimajoraxis;
Inclination(i) = x.catalog_TLEs{i}.inclination;
end

% 
% a = linspace(Earth.r,Earth.r+3000,1000);
% e = linspace(0,0.3,1000);
% i = linspace(90,180,1000);

J2 = 1.087*10^-3;%Equitorial Buldge

OmegaDot = 0.9856* (pi/180) * (1/24)*(1/60)*(1/60); %deg/day ->radians per second

f = @(a,e,i) (-3/2*J2*sqrt(mu)*Earth.r.^2) * cos(i)./(a.^(7/2)*(1-e.^2).^2) - OmegaDot; %deg/day BIG OMEGA DOT
interval = [Earth.r, Earth.r + 3000, 0, 0.3, 90, 180];
fimplicit3(f, interval)



f1 = @(x,y,z) x.^2+y.^2+z.^2;
fimplicit3(f, interval)
xlabel('a Semi Major(km)')
ylabel('e Eccentricity')
zlabel('i Inclination')
hold on
plot3(Semimajor,Eccentricity,Inclination,'.')
xlim([Earth.r Earth.r + 3000])

