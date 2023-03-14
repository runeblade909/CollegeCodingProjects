#PARETO CODE

#We will have a section for each high level objective (HLO), that results in a "score" or value that represents: Power Generated, Fuel Generated, ROI
#Iterate your trades and we'll combine all of the vectors and give each point scores on the three HLO's

#Libraries------------------------------------------------------------------------------------------
from csltk import cr3bp
from csltk.utilities import Earth, Moon, System
from csltk.correctors import SingleShooter

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
from math import log,exp,pow
from cmath import pi
import csv

import pandas as pd
import statistics
import numpy as np
import sys

import ROIPareto
import ROICalc
#import functions from other scripts in GroundTeamCode folder
sys.path.append('GroundTeamCode/')
import ExcavationCode
import thermal_pareto
import Process
import TransportCode

#import functions from other scripts in power_preliminary_models folder
sys.path.append('power_preliminary_models/')
import FullPowerSys

#Useful variables-----------------------------------------------------------------------------------
cost_of_kg = 20 #Highly dependent on LV (this is based on starship values)
sp_cost = 500/0.0032; #cost per square meter of solar panel
        
##FUEL TRADES ------------------------------------------------------------------------------


## ------------------------------------------------------------------------
## This code is to evaluate:

#Excavation Rate
#Power Consumption
#Simplicity

#Code by Zak Reichenbach


print("----------------------------------------------------------------")
print()
print()

ExcavationAnalysis = ExcavationCode.ExcavationAnalysis()
ExcavationAnalysis.excavator_productivity()

print(ExcavationAnalysis.Rho,"Density\n")
print(ExcavationAnalysis.PowerBucket,"Watt Power to operate Bucket Loader")
print(ExcavationAnalysis.ExcMassBucket, "Kg Bucket Material Mass")
print(ExcavationAnalysis.ExcRateBucket, "Kg/hr Bucket Loader Excavation Rate\n")


print(ExcavationAnalysis.PowerDrum,"Watt Power to operate Drum")
#print(ScoopSize, "Drum Scoop Size")
print(ExcavationAnalysis.ExcMassDrum, "Kg Drum Material Mass")
print(ExcavationAnalysis.ExcRateDrum, "Kg/hr Drum Excavation Rate\n")


print(ExcavationAnalysis.PowerConv, "Watts Power to operate Belt Conveyor")
print(ExcavationAnalysis.MassBin, "Kg Mass Holdable in storage")
print(ExcavationAnalysis.ExcMassConv, "Kg Mass Bucket Conveyor Material Mass")
print(ExcavationAnalysis.ExcRateConv, "Kg/hr Bucket Conveyor Rate\n")

print("\n\nNow for Excavator Productivity\n")

print("Bucket Productivity = ",ExcavationAnalysis.QBucket," m^3/hr\nDrum Productivity = ",ExcavationAnalysis.QDrum," m^3/hr\nConveyor Productivity = ",ExcavationAnalysis.QConv," m^3/hr\n")

print("Bucket Kg/hr = ",ExcavationAnalysis.RateBucket, "Drum Kg/hr = ",ExcavationAnalysis.RateDrum,"Conv Kg/hr = ",ExcavationAnalysis.RateConv)

print()

print("Final Vectors:\n",ExcavationAnalysis.BucketVec,"\n",ExcavationAnalysis.DrumVec,"\n",ExcavationAnalysis.ConvVec,"\n")

print("----------------------------------------------------------------")


## -------------------------------------------------------------------------


#Processing functions -----------------------------------------------------------------


print("\n\nNow for the Processing model for the Bucket, Drum, and Conveyor Systems\n\nBucket-------------------------------\n")

#Metals,Power,Cost,Fuel,Time,Tt1000

scale = [0]*3 #for scaling initial storage costs
ProcCost = [0]*3 #[3][2] array with all processing costs
ProcTime = [0]*3 #[3][2] array with all processing times

ProcBuckMetal,ProcBuckPower,ProcCost[0],ProcBuckFuel,ProcTime[0],ProcBuckTime1000,scale[0] = Process.Processing(ExcavationAnalysis.BucketVec[0])

ProcMethBucket = np.array([ProcBuckFuel[0],ProcBuckPower[0],0])
ProcHydrBucket = np.array([ProcBuckFuel[1],ProcBuckPower[1],0])
print("Drum-------------------------------\n")

ProcDrumMetal,ProcDrumPower,ProcCost[2],ProcDrumFuel,ProcTime[2],ProcDrumTime1000,scale[2] = Process.Processing(ExcavationAnalysis.DrumVec[0])

ProcMethDrum = np.array([ProcDrumFuel[0],ProcDrumPower[0],0])
ProcHydrDrum = np.array([ProcDrumFuel[1],ProcDrumPower[1],0])
print("Conv-------------------------------\n")

ProcConvMetal,ProcConvPower,ProcCost[1],ProcConvFuel,ProcTime[1],ProcConvTime1000,scale[1] = Process.Processing(ExcavationAnalysis.ConvVec[0])

ProcMethConv = np.array([ProcConvFuel[0],ProcConvPower[0],0])
ProcHydrConv = np.array([ProcConvFuel[1],ProcConvPower[1],0])

#print(BucketVec[0],"\n",DrumVec[0],"\n",ConvVec[0],"\n")

#print(ProcBuckMetal)


print("Final Processing Vectors\n\n",ProcMethBucket,"\n",ProcHydrBucket,"\n",ProcMethDrum,"\n",ProcHydrDrum,"\n",ProcMethConv,"\n",ProcHydrConv,"\n")
print("----------------------------------------------------------------")


#print(ProcBuckPower,"\n",ProcDrumPower,"\n",ProcConvPower,"\n")

#Thermal and Material Trades -----------------------------------------------------------------------------------------------

# Description: Thermal and Material Pareto
# Usage: python thermal_pareto.py
# Input: Thermal_Hardware.csv
# Output: thermal_pareto.png

#efficiency = mat out / power consumption
# for thermal this is thermal usage over power cosumption

# Assumptions: Similar rover to HAKUTO's SORATO rover
# https://ttu-ir.tdl.org/bitstream/handle/2346/86273/ICES-2020-209.pdf?sequence=1
# https://ocw.mit.edu/courses/16-851-satellite-engineering-fall-2003/e3a84cc153960fff8d55480fe228bbcc_l23thermalcontro.pdf

# rover will deal with Qrad, Qalbedo, Qsolar, Pdiss, Qir_moon
#    Qrad
#    Qalbedo
#    Qsolar
#    Qir_moon
material_data = thermal_pareto.read_excel('GroundTeamCode/Thermal_Hardware.csv')
Thermal_Analysis = thermal_pareto.ThermalAnalysis()
Thermal_Analysis.add_materials(material_data)
for i in range(len(Thermal_Analysis.materials[0])):
    Thermal_Analysis.solve_efficiency(i)
Thermal_Analysis.ranking = sorted(enumerate(Thermal_Analysis.efficiency), key = lambda i: i[1])

print("\nThermal Final Output Vector:\n\n")

#print(price)
#print("\n",eff)

Ranking = np.zeros((i,3))

#print(Ranking)
length = i
i = 0
#RoiScore = Profit
#PowerScore = a b c d e f
#FuelScore = of the -6 (Bucket, Drum, Conv) (Methane, Hydrogen)

#[Fuel, Power, ROI,ROIScore,PowerScore,FuelScore]


for i in range(0,length):
    Ranking[i,1] = Thermal_Analysis.efficiency[i]
    Ranking[i,2] = Thermal_Analysis.cost[i]
    Ranking[i,0] = 0

#print(Thermal_Analysis.ranking,"\n\n\n",price)

print(Ranking)

## Transport Code  ---------------------------------------------------

## Andrew Sapuppo 





# Transportation Code 
#_________________________________________________________________________________________________________________ 

#_________________________________________________________________________________________________________________ 

## Types of product delivery studies   

# DELIVERABLES: 

# (1) Extra cost for The Ice Box per sale [$] 

# (2) Extra ONE TIME cost for The Ice Box [$] 

# (3) Extra cost for customer per sale [$] 

# (4) Extra ONE TIME cost for customer [$] 

# (5) Amount of product that can be transported per sale [kg] 

# (6) If method of getting product to orbit is reusbale how much [kg]'s of fuel 

# or [J]'s of energy needed per sale 


# There are three types of transporation for product 

# (1) Customer picks up product on surface of moon at processing site 

# (2) We deliver to the customer in orbit 

#_________________________________________________________________________________________________________________ 

#_________________________________________________________________________________________________________________ 

## Assumptions 

# (1) Assumes all costs are equal to price to bring kg to lunar surface. Costs does not 

# account for bulding/developing vehicles or infrastructure OR include 

# buying a vehicle like Starship 

# (2) Since most of spin launch's infrastructure is for vacuum management assumes all other nessecary infrastructure is 100 metric tons 

# (3) To account for Moon's gravity assume we can multiply the amount of 

# mass spin lauch can lauch into orbit on Earth by this ratio to get moon equilvalent g/gm 

# (4) The extra mass needed to control the package of product after 

# launched to orbit from spin launch is assumed to be 200 kg since spin launch 

# currently offers this as the max weight satellite it can launch so assume 

# our package uses all that mass for trajectory/attitude control  

# (5) Assumes Starship can only take its minimum expected payload capacity 

# to the lunar surface. Its payload capacity is 100-150 metric tons 

# (6) Assumes Starship's propellant capacity is fully used to get to the 

# lunar surface 

# (7) To account for Moon's gravity assume we can divide the amount of 

# Starship propellant needed to leave Earth by this ratio to get moon equilvalent g/gm 

#_________________________________________________________________________________________________________________ 

#_________________________________________________________________________________________________________________ 


## TO DO NEED TO ADD OTHER RATES FOR FUEL, POWER, And excavation method

print("\n\n Transport Final Output Vectors\n\n")

TransportAnalysis = TransportCode.TransportAnalysis()


Time1000s = [ProcBuckTime1000[0],ProcDrumTime1000[0],ProcConvTime1000[0]]
Fuels = [ProcMethBucket[0],ProcHydrBucket[0],ProcMethDrum[0],ProcHydrDrum[0],ProcMethConv[0],ProcHydrConv[0]]

ParetoVec = np.zeros((3,3))

ParetoVector1 = []
ParetoVector2 = []
ParetoVector3 = []
ParetoVector4 = []
ParetoVector5 = []


print(Time1000s,"\n\n\n",Fuels,"\n\n\n")

i = 0

#ParetoVector1 = TransportAnalysis.ParetoVector1(Time1000s[0],Fuels[0])
#ParetoVector11 = TransportAnalysis.ParetoVector1(Time1000s[1],Fuels[1])

#print(Time1000s[1],"\n\n",Fuels[1],"\n\n",ParetoVector1,"\n\n")


#i loop is for each other the two methods (time for hydrogen and carbothermal), the j loop the different 

i = 0
counter = 0
for j in range(0,6):

        print("\n",i,j,"\n")
        ParetoVector1.append(TransportAnalysis.ParetoVector1(Time1000s[i],Fuels[j]))
        ParetoVector2.append(TransportAnalysis.ParetoVector2(Time1000s[i],Fuels[j]))
        ParetoVector3.append(TransportAnalysis.ParetoVector3(Time1000s[i],Fuels[j]))
        ParetoVector4.append(TransportAnalysis.ParetoVector4(Time1000s[i],Fuels[j]))
        ParetoVector5.append(TransportAnalysis.ParetoVector5(Time1000s[i],Fuels[j]))
        counter = counter + 1
        if counter == 2:
            counter = 0
            i = i+1
        




print(ParetoVector1," \n\n",ParetoVector2," \n\n",ParetoVector3," \n\n",ParetoVector4," \n\n",ParetoVector5)



#---------------------------------------------------------------------------


##POWER TRADES--------------------------------------------------------------------------------------

#Full power system function
#Contains all of the trades that go into the power system, and output the power generated, fuel generated, and cost and lifetime
#Inputs: 
#Dist_Pwr = power required from the distribution system
#Dist_energy = energy required from the distribution system
#Processing_Pwr = power required from the processing system
#Processing_Pwr = power required from the processing system
#a = battery trade variable (2 options)
#b = power hybrid architecture variable (4 options)
#c = RTG trade (2 option) hard set to 1
#d = Orbit trade (4 options) Might need to vary same ROI score
#e = Orbit Battery 0-1
#f = Recieving antenna diameter 0-2


#Dist_Pwr = power required from the distribution system
#Dist_energy = energy required from the distribution system
#Processing_Pwr = power required from the processing system
#Processing_Pwr = power required from the processing system
#a = battery trade variable (2 options)
#b = power hybrid architecture variable (4 options) 0- all orbital 3- all ground
#c = RTG trade (2 option) - num of rtgs
#d = Orbit trade (4 options) - different orbits
#e = Orbit Battery 0-1 - Orbital battery type
#f = Recieving antenna diameter 0-2 - 20-25-30


#Assuming we are always using a spinlaunch
#print("\n\n",ParetoVector5[5],"\n\n")

x = ParetoVector5[5]

Dist_Pwr = x[2]/(12*3600)
Dist_energy = x[2]

a = [0,1]
b = [0,1,2,3]
c = [0,1]
d = [0,1,2,3]
e = [0,1]
f = [0,1,2]
q = [0,1,2,3,4,5]

#Power,Fuel,Cost,Efficiency,Lifetime

loops = len(a)*len(b)*len(c)*len(d)*len(e)*len(f)*len(q)
length = 2304

Power = np.zeros((length,1))
Fuel = np.zeros((length,1))
Cost = np.zeros((length,1))
Efficiency = np.zeros((length,1))
Lifetime = np.zeros((length,1))
PowerCode = []


ProcessingPower = [ProcBuckPower[0],ProcBuckPower[1],ProcConvPower[0],ProcConvPower[1],ProcDrumPower[0],ProcDrumPower[1]]
ProcessingEnergy = ProcessingPower * (12*3600)
header = ['Power','Fuel','Cost','Efficiency','Lifetime','Power Code']


## Loop through power model and comine all options

iterations = 0

counter = 0

with open('powerLESS.csv', 'w') as file:
    writer = csv.writer(file)
    writer.writerow(header)

    for q in range(0,1): #Power 1-6 processing for excavation
        Processing_Pwr = ProcessingPower[q]
        Processing_energy = ProcessingEnergy[q]
        for i in range(0,1): #Battery 1-2
            a1 = a[i]
            #print(i)
            for j in range(0,1): # Hybrid Architecture ONLY ONE THAT CHANGES INFRASTRUCTURE COST 1-4
                b1 = b[j]
                for k in range(0,1): # RTG's 1-2
                    c1 = c[k]
                    for z in range(0,1): #Orbits 1-4
                        d1 = d[z]
                        for y in range(0,1): # Orbit Battery 1-2
                            e1 = e[y]
                            for x in range(0,1): # Reciever Diameter 1-3
                                f1 = f[x]
                                print("\nThis is loop number:",counter+1,"/",loops,"\n Values for loop PowerModel:",Processing_Pwr," a:",a1," b:",b1," c:",c1," d:",d1," e:",e1," f:",f1,"\n\n")
                                [Power[counter],Fuel[counter],Cost[counter],Efficiency[counter],Lifetime[counter]] = FullPowerSys.Power_Trades(Dist_Pwr, Dist_energy, Processing_Pwr, Processing_energy, a1, b1, c1, d1, e1, f1)
                                #print(Power,Fuel,Cost,Efficiency,Lifetime)
                                Code = str((q,a1,b1,c1,d1,e1,f1))
                                PowerCode.append(Code)
                                #print(float(Power[counter]),float(Fuel[counter]),float(Cost[counter]),float(Efficiency[counter]),float(Lifetime[counter]))
                                writer.writerow([Power[counter],Fuel[counter],Cost[counter],Efficiency[counter],Lifetime[counter],Code])
                                counter = counter + 1
                                #print("\n",i," ",j," ",k," ",z," ",y," ",x,"\n")

                            #if counter == 384:
                            #    iterations+=counter
                            #    writer.writerow(["-","-","-","-","-","-"])
                            #else:
                            #    iterations += counter


                            x = x+1
                        y = y+1
                    z = z+1
                k = k+1
            j = j+1
        i = i+1
    q = q + 1

file.close()


##https://www.geeksforgeeks.org/python-pandas-dataframe/

# Read in complete data set
data = pd.read_csv('power.csv')

print(ProcMethBucket,"\n",ProcHydrBucket,"\n",ProcMethDrum,"\n",ProcHydrDrum,"\n",ProcMethConv,"\n",ProcHydrConv,"\n")
print(data)
#print(data[["Power Code"]])
print(data["Power"][:]) 
#print("\n","\n",Power,"\n",Fuel,"\n",Cost,"\n",Efficiency,"\n",Lifetime,"\n")

OverAllVector = np.zeros((length,3))



PowerVec = np.zeros((length,1))

for k in range(0,len(Power)):
    PowerVec[k] = [Power[k]]




k=0

for i in range(0,length):
    #384
    #print("\n",k,"\n")

    ### This is a catch for the last index ####
    if k ==  6:
        k = 5


    OverAllVector[i,0] = Fuels[k]
    counter = counter + 1

    if counter == 384:
        counter = 0
        k = k+1

for k in range(0,length):
    #print(data["Power"][k])
    
   # counter = counter + 1
   # if counter == 384:
   #     counter = 0
   #     k = k+1

    if data["Power"][k][1:-1] == "":
        k = k+1

    OverAllVector[k,1] = float(data["Power"][k][1:-1])

        
#print("\n\n",OverAllVector,"\n\n")



#print("Final Power Vector:\n",PowerVec)


##ROI TRADES -----------------------------------------------------------------------------

# Pareto Code
# By: Amir Tillis
# Last Updated: 11/7/22
# Description: This code is intended to find the 2D pareto surface of
# revenue and infastructure 

#Code finds ROI score for 384 power, 5 tranport solutions, and 6 combined excavation/processing solutions


## Now we need to make all the datapoints from the 2304 choices go through the ROI Model




Profit = ROICalc.ROI(data,ProcMethBucket[0],ProcHydrBucket[0],ProcMethDrum[0],ProcHydrDrum[0],ProcMethConv[0],ProcHydrConv[0],ProcBuckMetal,ProcDrumMetal,ProcConvMetal,scale,ProcCost,ProcTime)

#print("\n\n",len(Profit),"\n",len(ROICode),"\n",len(Revenue_vs_Cost),"\n",len(Infrastructure_Cost),"\n")

ROIVector = Profit

#1-5 change in transport
#1-20 changes in power
#1-120 changes in ground method

##print(ROIVector)

#print(Infrastructure_Cost)
#print(Profit)

#Combining for ultimate results for plotting! --------------------------------

##### REMEMBER PARETO VEC FROM TRANSPORT MODEL ##############

ProcessingExcavation = np.zeros((6,3))

ProcessingExcavation[0,:] = ProcMethBucket
ProcessingExcavation[1,:] = ProcHydrBucket
ProcessingExcavation[2,:] = ProcMethDrum
ProcessingExcavation[3,:] = ProcHydrDrum
ProcessingExcavation[4,:] = ProcMethConv
ProcessingExcavation[5,:] = ProcHydrConv



#Power, Fuel, ROI
#print(ProcessingExcavation, "\n")
#ROI, Mass per hr, Power Per Sale
#print(ParetoVector1,"\n",ParetoVector2,"\n",ParetoVector3,"\n",ParetoVector4,"\n",ParetoVector5,"\n")
#Power, Fuel, ROI
#print(PowerVec,"\n")
#all 120 profit options, outline below
#print(Profit,"\n")
#Power code every different power option
#print(PowerCode,"\n")

#print(ROICode,"\n")

# for ROI loop
    #Only 6 Ground Options
    #Each Individual Ground Option has 7 Revenue Options (Wrapped into one total cost value), and 12 infastructure
    #cost options. Need to Compare each Revenue Option to 20 Infastrcture Cost (4 power options 5 transport)
    #1-5 change in transport
    #1-20 changes in power
    #1-120 changes in ground method

#overall_count = 0

#FinalVec = np.zeros((120,3))

#for i in range(0,6): #Ground
#        counter=0
#        #print("i:",i)
#        for j in range(0,4): #Power
#            #print("j:",j)
#            for h in range(0,5): #Transportation
#                
#                if k == 2: # if even PROCESSING METHOD
#                    k = 0
#
#                FinalVec[overall_count,:] = [ProcessingExcavation[i,0]+PowerVec[j,0],ProcessingExcavation[i,1]+TransportVector[h,1],Profit[overall_count]] #,PowerCode[j],ROICode[overall_count]]
#
#                counter += 1 
#                overall_count += 1
#    
#        k += 1

#Power,Fuel,Profit

#print("\n","\n",FinalVec,"\n","\n")



#print("\n\n",Profit,"\n\n")


OverAllVector[:,2] = Profit[0:-1,0]

#print("\n\n",OverAllVector,"\n\n")


FinalVec = OverAllVector

TrimmedVec = np.zeros((1500,3))

counter = 0

for i in range(0,length):

    if OverAllVector[i,2] > 0:
        TrimmedVec[counter,:] = [OverAllVector[i,0],OverAllVector[i,1],OverAllVector[i,2]]
        counter = counter + 1


print("\n\n",TrimmedVec,"\n\n",counter,"\n\n")

Index = np.zeros((500,1))
MaxFuel = np.zeros((500,3))

counter = 0
for i in range(0,length):

     if OverAllVector[i,0] > 1400:

        Index[counter] = i
        MaxFuel[counter] = OverAllVector[i,:]
        counter = counter+1

print("\n\n",MaxFuel,"\n\n")



w = np.max(MaxFuel[:,2])

counter = 0

thespot = np.zeros((10,1))

for k in range(0,500):
    if w == MaxFuel[k,2]:
        thespot[counter] = Index[k]
        print(thespot[counter])
        counter = counter + 1



print("\n\n",w,"\n\n",thespot,"\n\n")


for k in range(0,10):
    print(data["Power Code"][thespot[k]])


#---------------------------------------------------------------------------------------------
#power initial costs
#GS_cost = A_ground*sp_cost + m_ground*cost_of_kg #Cost of ground solar sytsem
#print(GS_cost)
#OS_cost = m_orbital*cost_of_kg + num_of_sat*cost_of_sat#Cost of orbital solar system
#RTG_cost = 100000000*np.array([0,1,2,3])#Cost of RTG's


#Plot Pareto from the three vectors above
# Fixing random state for reproducibility


#def randrange(n, vmin, vmax):
#    """
#    Helper function to make an array of random numbers having shape (n, )
#    with each number distributed Uniform(vmin, vmax).
#    """
#    return (vmax - vmin)*np.random.rand(n) + vmin

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

n = 120

# For each set of style and range settings, plot n random points in the box
# defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
xs = FinalVec[:,1]
ys = FinalVec[:,0]
zs = FinalVec[:,2]
ax.scatter(xs, ys, zs)

ax.set_xlabel('Power [kW]')
ax.set_ylabel('Fuel [kg/hr]')
ax.set_zlabel('Profit [$]')

plt.show()

fig2 = plt.figure()
ax2 = fig2.add_subplot(projection='3d')

xs2 = TrimmedVec[:,1]
ys2 = TrimmedVec[:,0]
zs2 = TrimmedVec[:,2]
ax2.scatter(xs2,ys2,zs2)

ax2.set_xlabel('Power [kW]')
ax2.set_ylabel('Fuel [kg/hr]')
ax2.set_zlabel('Profit [$]')

plt.show()



print(len(OverAllVector),"\n",len(TrimmedVec))