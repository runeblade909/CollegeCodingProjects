# Pareto Code
# By: Amir Tillis
# Last Updated: 11/7/22
# Description: This code is intended to find the 2D pareto surface of
# revenue and infastructure 

import numpy as np
import sys

#
#clc
#clear
#

## Main 

## Power

#Orbital Solar Array

OS_LS=15;  #Years
OS_Eff=100; #W
OS_Cost=110000000; #USD

OS_USD_per_kJ=OS_Cost*(1/OS_Eff)*(1/OS_LS)*(1/31536000); #$/kJ


#Ground Solar Array


GS_LS=14;  #Years
GS_Eff=100; #kW
GS_Cost=90000000; #USD

GS_USD_per_kJ=GS_Cost*(1/GS_Eff)*(1/GS_LS)*(1/31536000); #$/kJ


#50/50 Orbital Ground Mix

OG_LS=14;  #Years
OG_Eff=100; #kW
OG_Cost=100000000; #USD

OG_USD_per_kJ=OG_Cost*(1/OG_Eff)*(1/OG_LS)*(1/31536000); #$/kJ

## Ground 

#Bucket Loader
BL_LS=10;  #Years
BL_Eff=11812.5; #kg/hr
BL_Cost=10000000; #USD

BL_USD_per_kg=BL_Cost*(1/BL_Eff)*(1/BL_LS)*(1/870); #$/kg

#Bucket Conveyor 

BC_LS=5;  #Years
BC_Eff=22950; #kg/hr
BC_Cost=40000000; #USD

BC_USD_per_kg=BC_Cost*(1/BC_Eff)*(1/BC_LS)*(1/870); #$/kg

#Drum Loader

DL_LS=8;  #Years
DL_Eff=36000; #kg/hr
DL_Cost=30000000; #USD

DL_USD_per_kg=DL_Cost*(1/DL_Eff)*(1/DL_LS)*(1/870); #$/kg

## Processing Technology

#Carbothermal Redux

CR_LS=5;  #Years
CR_Eff=121.9; #kg/hr
CR_Cost=6887460; #USD

CR_USD_per_kg=CR_Cost*(1/CR_Eff)*(1/CR_LS)*(1/870); #$/kg


# Hydrogen Redux

HR_LS=5;  #Years
HR_Eff=87.64; #kg/hr
HR_Cost=1141181; #USD

HR_USD_per_kg=HR_Cost*(1/HR_Eff)*(1/HR_LS)*(1/870); #$/kg

# Transportation

#Technology 1

T1_LS=15;  #Years
T1_Eff=83.33; #kg/hr
T1_Cost=2000000; #USD

T1_USD_per_kg=T1_Cost*(1/T1_Eff)*(1/T1_LS)*(1/870); #$/kg

#Technology 2

T2_LS=10;  #Years
T2_Eff=27.89; #kg/hr
T2_Cost=2000000; #USD

T2_USD_per_kg=T2_Cost*(1/T2_Eff)*(1/T2_LS)*(1/870); #$/kg

#Technology 3

T3_LS=5;  #Years
T3_Eff=35.19; #kg/hr
T3_Cost=294835.04; #USD

T3_USD_per_kg=T3_Cost*(1/T3_Eff)*(1/T3_LS)*(1/870); #$/kg

#Technology 4

T4_LS=5;  #Years
T4_Eff=35.19; #kg/hr
T4_Cost=28520876.19; #USD

T4_USD_per_kg=T4_Cost*(1/T4_Eff)*(1/T4_LS)*(1/870); #$/kg

#Technology 5

T5_LS=5;  #Years
T5_Eff=35.19; #kg/hr
T5_Cost=294835040.8; #USD

T5_USD_per_kg=T5_Cost*(1/T5_Eff)*(1/T5_LS)*(1/870); #$/kg






## Products

#Regolith
Regolith_Abd=1;
Regolith_Cost=25; #USD/kg

#Oxygen
Oxygen_Abd=.45;
Oxygen_Cost=11.8; #USD/kg

#Silicon
Silicon_Abd=.22;
Silicon_Cost=9.60; #USD/kg

#Iron
Iron_Abd=.13;
Iron_Cost=.08; #USD/kg

#Calcium
Calcium_Abd=.08;
Calcium_Cost=4.9; #USD/kg

#Aluminum
Aluminum_Abd=.07;
Aluminum_Cost=2.25; #USD/kg

#Magnesium
Magnesium_Abd=.05;
Magnesium_Cost=3.31; #USD/kg


## Revenue and COst Calculations

Ground_Eff=[BL_Eff, BC_Eff, DL_Eff] 
Ground_LS=[BL_LS,BC_LS,DL_LS]
Product_Abd=[Regolith_Abd,Oxygen_Abd,Silicon_Abd,Iron_Abd,Calcium_Abd,Aluminum_Abd,Magnesium_Abd]
Product_Cost=[Regolith_Cost,Oxygen_Cost,Silicon_Cost,Iron_Cost,Calcium_Cost,Aluminum_Cost,Magnesium_Cost]

print("Ground_Eff")
print(Ground_Eff)
print("Ground_LS")
print(Ground_LS)
print("Product_Abd")
print(Product_Abd)
print("Product_Cost")
print(Product_Cost)


Revenue = np.zeros((3,7))


#print(float(Ground_Eff[0])*float(Ground_LS[0])*(870)*float(Product_Abd[3])*float(Product_Cost[3]))

count = 1

for i in range(0,3): #Ground
    for j in range(0,7): #Products
        #print(count)
        #count = count + 1
        #print(float(Ground_Eff[i])*float(Ground_LS[i])*(870)*float(Product_Abd[j])*float(Product_Cost[j]))
        Revenue[i,j] = float(Ground_Eff[i])*float(Ground_LS[i])*(870)*float(Product_Abd[j])*float(Product_Cost[j]);

print("Revenue")
print(Revenue)

Power_Eff=[OS_Eff,GS_Eff,OG_Eff];
Power_LS=[OS_LS,GS_LS,OG_LS];

Power_Cost=[OS_USD_per_kJ,GS_USD_per_kJ,OG_USD_per_kJ];
Ground_Cost=[BL_USD_per_kg,BC_USD_per_kg,DL_USD_per_kg];
Processing_Cost=[CR_USD_per_kg,HR_USD_per_kg];
Transportation_Cost=[T1_USD_per_kg,T2_USD_per_kg,T3_USD_per_kg,T4_USD_per_kg,T5_USD_per_kg];

Infastructure_Cost= np.zeros((3,30))

#print(Power_Eff)
#print(Power_LS)
#print(Power_Cost)
#print(Ground_Cost)
#print(Processing_Cost)
#print(Transportation_Cost)


overall_count = 0


for i in range(0,3): #Ground
    counter=0
    for j in range(0,3): #Power
        for k in range(0,2): #Processing
            for h in range(0,5): #Transportation
                Infastructure_Cost[i,counter] = float(Ground_LS[i])*float(Ground_Eff[i])*870*(float(Ground_Cost[i])+float(Processing_Cost[k])+float(Transportation_Cost[h]))+float(Power_LS[j])*float(Power_Eff[j])*31536000*float(Power_Cost[j])
                counter += 1 
                overall_count += 1


print("Infastructure_Cost")
print(Infastructure_Cost)

#Only 3 Ground Options

#Each Individual Ground Option has 7 Revenue Options, and 12 infastructure
#cost options. Need to Compare each Revenue Option to 12 Infastrcture Cost
#Options 

Revenue_vs_Cost = np.zeros((3*7*12,2)) #[[0 for col in range(3*7*12)] for row in range(2)]
counter = 0
for i in range(0,3): #Mining Option
    for j in range(0,7): #Revenue 
        for k in range(0,12): #Cost
        
                Revenue_vs_Cost[counter,:]=[Revenue[i,j],Infastructure_Cost[i,k]]; 
                counter=counter+1;
            #figure(1)
            #plot(Revenue(i,j),Infastructure_Cost(i,k),'o');
            #hold on
            #xlabel('LifeTime Revenue USD')
            #ylabel('Infastructure Cost USD')
            #title('LifeTime Revenue vs Infastrcture Cost')
print("Revenue v Cost")       
print(Revenue_vs_Cost)

Profit= Revenue_vs_Cost[:,0] - Revenue_vs_Cost[:,1]
print("Profit")
print(Profit)