%% Proposed cost savings from providing material on the moon.

% Zak Reichenbach

clc
clear all 
close all

%% Price per Kg

%Based off of Falcon 9 Launch Characteristics.

LaunchCost = 67000000; %$ Falcon 9 Cost

GTOMass = 8300; %kg Mass to GTO

MarsMass = 4020; % kg Mass to Mars

LLOMass = (GTOMass + MarsMass)/2; % kg Guess Mass to Moon

BlueGhostCost = 16000000; %$

PPKGtoLLO = (LaunchCost+BlueGhostCost)/LLOMass; % $/kg

% Rover 
RoverCost = 40000000; %$

MiningSpeed = 40; %kg/hr up to 100 kg/hr

DailySpeed = MiningSpeed*24; %kg/hr

YearlySpeed = MiningSpeed*24*365*.94; %kg/hr

CurrentLifetime = 2.5; %years

%Regolith Content

Oxygen = 0.45;
Silicon = 0.22;
Iron = .13;
Calcium = 0.08;
Aluminum = 0.07;
Magnesium = 0.05;
Other = 0.1;

SpTitanium = 0.075; % in the maria


%% Our Finances

OurCost = RoverCost*3 + BlueGhostCost + LaunchCost;

%Products created from mission
ProductMass = YearlySpeed * 0.8; % if we produce everying with 80% efficiency from mass
%kgs produced
OxygenMass = ProductMass * Oxygen;
SiliconMass = ProductMass * Silicon;
IronMass = ProductMass * Iron;
CalciumMass = ProductMass * Calcium; 
AluminumMass = ProductMass * Aluminum;
MagnesiumMass = ProductMass * Magnesium;
SpTitaniumMass = ProductMass * SpTitanium;
OtherMass = ProductMass * Other;

OurPPKgMoon = OurCost/ProductMass;

%Material Price
OxDensity = 1141; % kg/m^3

MASS = OxDensity * 0.00378541; % mass of liquid oxygen in 1 gallon

OxPrice = 1.65 ;% $/gal --> $/kg 4.3192 kg

OxPrice = OxPrice/MASS; % $/kg
SiPrice = (20000/7.19)/1000; % rmb/mt --> USD/kg
FePrice = 0.09735;% $/kg
TiPrice = 18; % $/kg

%Regular material Price for our production
OxDebit = OxPrice*OxygenMass; % $ value 
SiDebit = SiPrice*SiliconMass;% $
FeDebit = FePrice*IronMass; % $
TiDebit = TiPrice*SpTitaniumMass; % $

MaterialsSold = OxDebit + SiDebit + FeDebit + TiDebit; % Total money made selling product for earth prices

%Competitor Pricing
NumLaunch = ProductMass/LLOMass;

CompetitorCost = BlueGhostCost*ceil(NumLaunch) + LaunchCost*ceil(NumLaunch)+MaterialsSold*.7;

CompPPkgMoon = CompetitorCost/ProductMass;

% How much we would save people to launch what we produce in a year
Savings = CompetitorCost - OurCost;

%What we can charge (Our cost + Material Sold Margin)/
Sales = (OurCost + MaterialsSold*100)/ProductMass;


TotalRoverProduct = CompetitorCost * 2.5


%rover can generate X material 