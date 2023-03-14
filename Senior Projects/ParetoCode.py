#PARETO CODE

#We will have a section for each high level objective (HLO), that results in a "score" or value that represents: Power Generated, Fuel Generated, ROI
#Iterate your trades and we'll combine all of the vectors and give each point scores on the three HLO's

#Libraries------------------------------------------------------------------------------------------

from csltk import cr3bp
from csltk.utilities import Earth, Moon, System
from csltk.correctors import SingleShooter

import matplotlib.pyplot as plt
import math
import csv

import numpy
import pandas as pd
import statistics
import numpy as np
import matplotlib.cm as cm


#Useful variables-----------------------------------------------------------------------------------
cost_of_kg = 20 #Highly dependent on LV (this is based on starship values)
sp_cost = 500/0.0032; #cost per square meter of solar panel
        
##FUEL TRADES ------------------------------------------------------------------------------



## ------------------------------------------------------------------------
## This code is to evaluate

#Excavation Rate
#Power Consumption
#Simplicity

#Code by Zak Reichenbach

from cmath import pi

print()
print()


#Regolith Properties
Rho = 1.5 #g/cm^3
Rho = Rho *100**3/1000 #Kg.m^3

print(Rho,"Density\n")

#Power Constraints
MotorPower = 50 # Each motor takes 50W
TorqueMotor = 75 # Higher torque application motor

#Bucket loaded
VolBucket = .6*.5*.5  # Volume of bucket m^3

ExcMassBucket = Rho * VolBucket  #Bucket mass 

ExcRateBucket = ExcMassBucket *.75 * 60#Kg/hr

SimplicityBucket = 8

PowerBucket = TorqueMotor*2 #2 motor linear actuators to raise and move bucket 

print(PowerBucket,"Watt Power to operate Bucket Loader")
print(ExcMassBucket, "Kg Bucket Material Mass")
print(ExcRateBucket, "Kg/hr Bucket Loader Excavation Rate\n")

#Drum Excavator
r = .25 #m
L = 0.5 #m
VolDrum = pi*r**2 *L #m^3 drum volume

ExcMassDrum = Rho * VolDrum  #kg

ScoopSize = .1*.1*.2   #m^3 for 5 scoops across the length, 3 distinct rows (see NASA Rassor)

ExcRateDrum = 12 * (ScoopSize*15)*Rho*60 # Kg/hr

PowerDrum = TorqueMotor*4 #2 motors to move bucket and rotate bucket

print(PowerDrum,"Watt Power to operate Drum")
#print(ScoopSize, "Drum Scoop Size")
print(ExcMassDrum, "Kg Drum Material Mass")
print(ExcRateDrum, "Kg/hr Drum Excavation Rate\n")

#Bucket Conveyor
OneBucket = .5*.1*.1

VolBin = .7*.3*.5 #Volume of storage m^3 bin for conveyor belt regolith

MassBin = VolBin*Rho # Mass kg

VolConv = OneBucket * 10

ExcMassConv = Rho * VolConv/2

ExcRateConv = ExcMassConv * 12 * 60 #Kg/hr

PowerConv = MotorPower + 4*TorqueMotor # 1 linear actuator motor to lower to ground, 2 motors for rotating coveyor and positioning

print(PowerConv, "Watts Power to operate Belt Conveyor")
print(MassBin, "Kg Mass Holdable in storage")
print(ExcMassConv, "Kg Mass Bucket Conveyor Material Mass")
print(ExcRateConv, "Kg/hr Bucket Conveyor Rate\n")

print("\n\nNow for Excavator Productivity\n")

#I used this resource for productivity
#https://homesteady.com/13650496/how-to-calculate-excavator-productivity

#Q = (60*q*z*n*kf)/kl [cubic ft/hr]
#q = volume of bucket
#z = number of buckets
#n = RPM or speed of rotation of the rotor
#kf = filling factor of bucket (0-1)
#kl = soil-loosening factor ~1.2

QBucket = (60*VolBucket*1*1.5*1*.7)/1.2
QDrum = (60*ScoopSize*15*12*.85)/1.2
QConv = (60*OneBucket*10*12*.8)/1.2

print("Bucket Productivity = ",QBucket," m^3/hr\nDrum Productivity = ",QDrum," m^3/hr\nConveyor Productivity = ",QConv," m^3/hr\n")

RateBucket = QBucket * Rho
RateDrum = QDrum * Rho
RateConv = QConv * Rho

print("Bucket Kg/hr = ",RateBucket, "Drum Kg/hr = ",RateDrum,"Conv Kg/hr = ",RateConv)



## Final vector scores
BucketVec = np.array([RateBucket,PowerBucket,0])
DrumVec = np.array([RateDrum,PowerDrum,0])
ConvVec = np.array([RateConv,PowerConv,0])


print()

print("Final Vectors:\n",BucketVec,"\n",DrumVec,"\n",ConvVec,"\n")


## -------------------------------------------------------------------------


from math import log,exp,pow

#Processing functions -----------------------------------------------------------------


def MOredux(MO,y,x):
    #x = mol of O, y = mol of M
    M = y*MO
    CH4req = x*MO
    H2 = 2*x*MO
    CO = x*MO
    return CH4req,H2,CO,M

def H2redux(MO,y,x):
    #x = mol of O, y = mol of M
    M = y*MO
    H2req = x*MO
    H20 = x*MO
    return H2req,H20,M

def COredux(CO,H2):
    molC = 12.011
    molH = 1.008
    molO = 15.999
    opTemp1 = 1625+273.15
    opTemp2 = 425+273.15; #K <-- 698.15 K
    cpCO = 1.113*1000/2; #J/kg/K @700
    cpH2 = 14.6*1000; #J/kg/K @700
    
    mCO = CO*(molO+molC); #mass of CO
    mH2 = H2*2*molH; #mass of H2
    CH4 = H2/3;
    H2O = H2/3;
    COex = CO-H2O;
    
    #assuming 500 K loss from storage 1 to 2
    P = (mCO+mH2)*((mCO/(mCO+mH2))*cpCO+(mH2/(mCO+mH2))*cpH2)*(opTemp2-(opTemp1-500)) #pwr required to mantain optemp
    return CH4,COex,H2O,P

def electro(H2O):
    P = 237000 * H2O #J
    H2 = H2O
    O2 = H2O/2
    return H2,O2,P

# Processing Main

def Processing(regRate):

    # Mass of Regolith
    #Reg = 1000; # Regolith (kg)
    time = 1 #hr
    Reg = regRate*time
    m = Reg; # mass of regolith 
    
    Tday = 260; #lunar surface temp at equator during day
    Tnight = 140; #lunar surface temp at equator during night
    #time = Reg/(regRate*3600); #time since start of process
    
    Reg = regRate*time

    pPkg = 5000; #price per kg

    # Mass of Metaloxides (using average values)
    mSiO2 = 0.5*Reg; #      Silicon-dioxide (kg)
    mAl2O3 = 0.15*Reg; #    Aluminum-oxide (kg)
    mCaO = 0.1*Reg; #       Calcium-oxide (kg)
    mFeO = 0.1*Reg; #       Iron-oxide (kg)
    mMgO = 0.1*Reg; #       Magnisium-oxide (kg)
    mTiO2 = 0.05*Reg; #     Titanium-oxide (kg)

    # Molar mass of Metaloxides
    molH = 1.008;
    molC = 12.011;
    molSi = 28.086;
    molO = 15.999;
    molAl = 26.982;
    molCa = 40.078;
    molFe = 55.845;
    molMg = 24.305;
    molTi = 47.867;

    # Moles of Metaloxides
    SiO2 = 1/(molSi+2*molO)    *mSiO2; #   Silicon-dioxide (mol)
    Al2O3 = 1/(2*molAl+3*molO) *mAl2O3; #  Aluminum-oxide (mol)
    CaO = 1/(molCa+molO)       *mCaO; #    Calcium-oxide (mol)
    FeO = 1/(molFe+molO)       *mFeO; #    Iron-oxide (mol)
    MgO = 1/(molMg+molO)       *mMgO; #    Magnisium-oxide (mol)
    TiO2 = 1/(molTi+2*molO)    *mTiO2; #   Titanium-oxide (mol)


    #specific heats of regolith at 1370 and 1898.15 (J/kg/K)
    cpR1370 = exp(((3*pow(log(1370),3)) + (-54.4484*pow(log(1370),2)) + (306.8465*log(1370)) - 376.5795)/(pow(log(1370),2) + (-16.8077*log(1370)) + 87.3213));

    cpR1898 = exp(((3*pow(log(1898.15),3)) + (-54.4484*pow(log(1898.15),2)) + (306.8465*log(1898.15)) - 376.5795)/(pow(log(1898.15),2) + (-16.8077*log(1898.15)) + 87.3213));

    ## Carbothermal Reduction

    CH4 = [0]*6
    H2 = [0]*6
    CO = [0]*6

    CH4[0],H2[0],CO[0],Si = MOredux(SiO2,1,2);
    CH4[1],H2[1],CO[1],Al = MOredux(Al2O3,2,3);
    CH4[2],H2[2],CO[2],Ca = MOredux(CaO,1,1);
    CH4[3],H2[3],CO[3],Fe = MOredux(FeO,1,1);
    CH4[4],H2[4],CO[4],Mg = MOredux(MgO,1,1);
    CH4[5],H2[5],CO[5],Ti = MOredux(TiO2,1,2);

    CH4req = sum(CH4);
    H2tot = sum(H2);
    COtot = sum(CO);
    opTemp1 = 1625+273.15; #K
    k = pow(10,(-2126.1/(opTemp1+0.6439)));
    time = time + (0.7/k); 
    Ptot = m*cpR1898*(opTemp1-Tday); #pwr required to mantain optemp


    CH4,COtot,H2O,P = COredux(COtot,H2tot);
    H2tot = 0; #assuming all hydrogen converted to H20
    CH4req = CH4req - CH4;
    time = time + (5*3600); 
    Ptot = Ptot+P;

    H2e,O2tot,P = electro(H2O);
    H2tot = H2e;
    time = time + (3.2*3600); #assuming 3.2 hrs for electrolysis
    Ptot = Ptot+P;
    

    Power = [Ptot/time,0]; #J/s
    Time = [time/3600,0]; #hr
    Cost = [(CH4req*pPkg*(molC+4*molH))/Time[0],0]; #$
    Fuel = [2*(H2tot*2*molH) + O2tot*2*molO, 0]; #score

    print('\nMethane Reduction:\n\n');
    print(f'Total Power required:\t {Ptot} J \n');
    print(f'Total Hydrogen gas produced:\t {H2tot*2*molH} kg \n');
    print(f'Total Oxygen gas produced:\t {O2tot*2*molO} kg \n');
    print(f'Total Methane gas required:\t {CH4req*(molC+4*molH)} kg \n');
    print(f'Total Cost of process:\t $ {CH4req*pPkg*(molC+4*molH)} \n');
    print(f'Total Time for process:\t {time/3600} hr \n');
    print(f'Arbitrary rate:\t {Ptot/time} W \n');
    print('Total Metals produced: \n');
    print(f'\t Ti:\t {Ti*molTi} kg \n');
    print(f'\t Si:\t {Si*molSi} kg \n');
    print(f'\t Mg:\t {Mg*molMg} kg \n');
    print(f'\t Al:\t {Al*molAl} kg \n');
    print(f'\t Fe:\t {Fe*molFe} kg \n');
    print(f'\t Ca:\t {Ca*molCa} kg \n');

    ## Hydrogen Reduction

    H2 = [0]*6
    H2O = [0]*6

    H2[0],H2O[0],Si = H2redux(SiO2,1,2);
    H2[1],H2O[1],Al = H2redux(Al2O3,2,3);
    H2[2],H2O[2],Ca = H2redux(CaO,1,1);
    H2[3],H2O[3],Fe = H2redux(FeO,1,1);
    H2[4],H2O[4],Mg = H2redux(MgO,1,1);
    H2[5],H2O[5],Ti = H2redux(TiO2,1,2);

    H2req = sum(H2);
    H2Otot = sum(H2O);
    opTemp1 = 1370; #K
    k = pow(10,(-2126.1/(opTemp1+0.6439)));
    time = time + (0.7/k); 
    Ptot = m*cpR1370*(opTemp1-Tday); #pwr required to mantain optemp

    H2e,O2tot,P = electro(H2Otot);
    H2tot = H2e-H2req;
    time = time + (3.2*3600); #assuming 3.2 hrs for electrolysis
    Ptot = Ptot+P;

    Power[1] = Ptot/time; #J/s
    Time[1] = time/3600; #hr
    Cost[1] = (H2req*(2*molH)*pPkg)/Time[1]; #$
    Fuel[1] = 2*(H2tot*2*molH) + O2tot*2*molO; #score


    print('\nHydrogen Reduction:\n\n');
    print(f'Total Power required:\t {Ptot} J \n');
    print(f'Total Hydrogen gas produced:\t {H2tot*2*molH} kg \n');
    print(f'Total Oxygen gas produced:\t {O2tot*2*molO} kg \n');
    print(f'Total Hydrogen gas required:\t {H2req*(2*molH)} kg \n');
    print(f'Total Cost of process:\t $ {H2req*(2*molH)*pPkg} \n');
    print(f'Total Time for process:\t {time/3600} hr \n');
    print(f'Arbitrary rate:\t {Ptot/time} W \n');
    print('Total Metals produced: \n');
    print(f'\t Ti:\t {Ti*molTi} kg \n');
    print(f'\t Si:\t {Si*molSi} kg \n');
    print(f'\t Mg:\t {Mg*molMg} kg \n');
    print(f'\t Al:\t {Al*molAl} kg \n');
    print(f'\t Fe:\t {Fe*molFe} kg \n');
    print(f'\t Ca:\t {Ca*molCa} kg \n');
    
    return(Power,Cost,Fuel,Time)

print("\n\nNow for the Processing model for the Bucket, Drum, and Conveyor Systems\n\nBucket-------------------------------\n")

ProcBuckPower,ProcBuckCost,ProcBuckFuel,ProcBuckTime = Processing(BucketVec[0])
ProcMethBucket = np.array([ProcBuckFuel[0],ProcBuckPower[0],0])
ProcHydrBucket = np.array([ProcBuckFuel[1],ProcBuckPower[1],0])
print("Drum-------------------------------\n")

ProcDrumPower,ProcDrumCost,ProcDrumFuel,ProcDrumTime = Processing(DrumVec[0])
ProcMethDrum = np.array([ProcDrumFuel[0],ProcDrumPower[0],0])
ProcHydrDrum = np.array([ProcDrumFuel[1],ProcDrumPower[1],0])
print("Conv-------------------------------\n")

ProcConvPower,ProcConvCost,ProcConvFuel,ProcConvTime = Processing(ConvVec[0])
ProcMethConv = np.array([ProcConvFuel[0],ProcConvPower[0],0])
ProcHydrConv = np.array([ProcConvFuel[1],ProcConvPower[1],0])

#print(BucketVec[0],"\n",DrumVec[0],"\n",ConvVec[0],"\n")


print("Final Processing Vectors\n\n",ProcMethBucket,"\n",ProcHydrBucket,"\n",ProcMethDrum,"\n",ProcHydrDrum,"\n",ProcMethConv,"\n",ProcHydrConv,"\n")


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


# Analyze:
#    - Prevent heat escape -- gold paint, insulation, "aerogel"
#    - Keep rover warm -- heater
#    - heat rejection -- white paint


# Rover Design Temps
# Not exceed -40degC to +40degC (mars)

# rover will deal with Qrad, Qalbedo, Qsolar, Pdiss, Qir_moon
#    Qrad
#    Qalbedo
#    Qsolar
#    Qir_moon
# solar power flux: 1367 W/m2 +- 3.5%
# sink temp is 2.53 K



def read_excel(file_name):
	data = pd.read_csv(file_name)
	return data

# function that searches through list linearly and returns desired idx
def linear_search(L, k):  # for small array you won't look at more than once
    for i in range(len(L)):
        if k == L[i]:
            return i
    return -1

class ThermalAnalysis:

	def __init__(self):
		# active or passive
		#self.type = type
		self.Stefan_Bolz_Const = 5.670 * pow(10,-8) # W/m^2-K^4
		self.Temp_Sink = 2.53 # Kelvin
		self.Solar_Flux = 1418 #G_s [W/m^2]
		self.Lunar_Albedo = 0.18  #min: 0.076, max:0.297
		# assuming dimensions of 529 x 602 x 372  [mm]
		self.Rover_Surface_Area = []#1478.38 #[m]
		self.Rover_Altitude = 0.2 #[m]
		self.Energy_Flux_Lunar_Surface = 0.031 # +- 20% [W/m^2]
		self.Solar_Power_Flux = 1367 # +-3.5% [W/m^2]
		# https://www.ri.cmu.edu/pub_files/pub1/deans_matthew_1998_1/deans_matthew_1998_1.pdf
		self.Solar_Incidence_Angle = 1.5 #[deg]
		self.P_disp = 32.6 # [W] from ICES paper
		self.P_disp_safe = 11.2 # [W]

		self.materials = []
		self.solar_absorptivity = []
		self.emmissivity = []
		self.power_requirement = []
		# dollars per kg
		# per unit - https://www.edmundoptics.com/p/10mm-dia-x-2mm-lambda4-first-surface-mirror/1199/
		# 
		self.cost = []
		self.surface_temp_min = []
		self.surface_temp_max = []
		self.surface_temp_avg = []
		self.thermal_balance = []
		self.efficiency = []
		self.materials_string = []
		self.ranking = []
		
	# material must be a panda containing info about the material
	def add_materials(self, material_data):
		self.materials.append(material_data.Material)
		self.solar_absorptivity.append(material_data.Solar_Absorptivity)
		self.emmissivity.append(material_data.Emmissivity)
		self.power_requirement.append(material_data.Power_Req)
		for i in material_data.Cost:
			self.cost.append(float(i))
		#self.cost = self.cost.replace(',','')
		#self.cost.append([int(cost) for cost in material_data.Cost])

		self.surface_temp_min.append(material_data.Surface_Temp_Min + 273)
		self.surface_temp_max.append(material_data.Surface_Temp_Max + 273)
		for i in range(len(self.surface_temp_min[0])):
			self.surface_temp_avg.append(statistics.mean([self.surface_temp_min[0][i], self.surface_temp_max[0][i]]))
			self.materials_string.append(str(self.materials[0][i]))
		return 1

	#def get_materials(self):
	#	return selfm
	def solve_thermal_balance(self, idx):
		# modify SA to solve thermal balance
		flag = False
		SA = 100 # [m]
		while flag == False:
			Q_rad = self.emmissivity[0][idx] * self.Stefan_Bolz_Const * SA * pow(self.surface_temp_avg[idx] - self.Temp_Sink, 4)
			
			K_a = 0.664 + 0.521 * self.Rover_Altitude + 0.203 * pow(self.Rover_Altitude, 2) #from MIT slides
			Q_albedo = self.Solar_Flux * self.Lunar_Albedo * SA * self.solar_absorptivity[0][idx] * K_a * pow(math.sin(self.Rover_Altitude), 2)
			
			Q_IR_moon = self.Energy_Flux_Lunar_Surface * pow(math.sin(self.Rover_Altitude), 2) * SA * self.emmissivity[0][idx]

			Q_solar = self.Solar_Power_Flux * SA * self.solar_absorptivity[0][idx] * math.cos(self.Solar_Incidence_Angle)

			thermal_balance = (-1 * Q_rad) + Q_albedo + Q_IR_moon + Q_solar + self.P_disp
			#print(thermal_balance)
			if thermal_balance < 50 and thermal_balance > -50:
				flag = True
				self.thermal_balance.append(thermal_balance)
				self.Rover_Surface_Area.append(SA)
			SA -= 0.05

		return 1

	def solve_efficiency(self, idx):
		# assumption: using T_surroundings of Sink Temp ?
		
		ThermalAnalysis.solve_thermal_balance(self, idx)

		self.efficiency.append(abs(self.thermal_balance[idx] / self.power_requirement[0][idx]))

		return self.efficiency[idx]


Eff = np.zeros((1,13))

if __name__ == '__main__':
	material_data = read_excel('Thermal_Hardware.csv')

	Thermal_Analysis = ThermalAnalysis()
	Thermal_Analysis.add_materials(material_data)

	
	Efficiency = []
	Labels = []
	for i in range(len(Thermal_Analysis.materials[0])):
		Efficiency.append(Thermal_Analysis.solve_efficiency(i))
	Thermal_Analysis.ranking = sorted(enumerate(Efficiency), key = lambda i: i[1])
	for i in range(len(Thermal_Analysis.materials_string)):
		Labels.append(Thermal_Analysis.materials_string[Thermal_Analysis.ranking[i][0]])
	fig, ax = plt.subplots()
	fig.set_size_inches(10, 5)
	# TODO: GET QUOTES 
	# PRICE IS RANDOM ARRAY NOT REAL VALUES
	price = np.random.uniform(1, 10, size=len(Thermal_Analysis.materials[0]))

	colors = cm.rainbow(np.linspace(0, 1, len(Thermal_Analysis.thermal_balance)))
	scatter = ax.scatter(Thermal_Analysis.thermal_balance, Thermal_Analysis.power_requirement, c = list(zip(*Thermal_Analysis.ranking))[0], s = 0.3*(price*3)**2, cmap = 'Spectral')
	legend1 = ax.legend(*scatter.legend_elements(num=len(Thermal_Analysis.materials[0])), loc = 'upper left', title = "Efficiency", bbox_to_anchor=(1, 1))
	ax.add_artist(legend1)
	kw = dict(prop="sizes", num=5, color=scatter.cmap(0.7), fmt="$ {x:.2f}", func=lambda s: np.sqrt(s/.3)/3)
	legend2 = ax.legend(*scatter.legend_elements(**kw), loc="upper right", title="Price")
	ax.set_xlabel("Thermal Balance")
	ax.set_ylabel("Power Required")
	ax.set_title("Thermal Efficiency")

	
	plt.savefig('thermal_pareto')
	print("Materials Ranked by Efficiency: ")
	print("\n")
	i = 0 
	for eff in Thermal_Analysis.ranking:
		print(str(i+1)+ " : " + Labels[i] + " Efficiency: " + str(eff[1]) + ", Thermal Balance: " + str(Thermal_Analysis.thermal_balance[i]) + ", Max SA: " + str(Thermal_Analysis.Rover_Surface_Area[i]))
		i += 1

	i = 0 
	print("\n")
	print("Price (Currently Rand Values) : ")
	for i in range(len(Thermal_Analysis.ranking)):
		print(str(i+1)+ " : " + Labels[i] + " Price: $" + str(price[i]))


print("\nFinal Output vector\n\n")


#print(price)
#print("\n",eff)


Ranking = np.zeros((i,3))

#print(Ranking)
length = i
i = 0






for i in range(0,length):
    Ranking[i,1] = Thermal_Analysis.efficiency[i]
    Ranking[i,2] = price[i]
    Ranking[i,0] = 0

#print(Thermal_Analysis.ranking,"\n\n\n",price)

print(Ranking)

## Transport Code  ---------------------------------------------------

print("\n\n Transport Final Output Vectors\n\n")

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

## Constants  

g = 9.81                   # Earth Gravitational constant [m/s^2] 

gm = 1.625                 # Lunar Gravitational constant [m/s^2] 

lb2kg = 2.20462262         # Divide lb by this to get [kg] 

mph2ms = 2.23693629        # Divide mph by this to get [m/s] 


#_________________________________________________________________________________________________________________ 

#_________________________________________________________________________________________________________________ 

## Saturn V takes Apollo 11/Blue Origin type lander to Moon (Assumes only ascent vehicle is reusable) 
 

p = 20000                                  # Price to launch one kg to the moon [$] (BASED on Apollo 11) 


# Apollo 11/Blue Origin type lander 

m11w = 32500 / lb2kg             # Wet mass of apollo 11 lunar lander (w/ propellant and crew) [kg] 

m11d = 9000 / lb2kg              # Dry mass of apollo 11 lunar lander (w/o propellant and crew) [kg] 
 

mA11w = 4700                     # Wet mass of apollo 11 lunar ascent vechicle (w/ propellant and crew) [kg] 

mA11d = 2150                     # Dry mass of apollo 11 lunar ascent vechicle (w/o propellant and crew) [kg] 
 

mD11w = m11w - mA11w             # Wet mass of apollo 11 lunar descent vechicle (w/ propellant and crew) [kg] 

mD11d = m11d - mA11d             # Dry mass of apollo 11 lunar descent vechicle (w/o propellant and crew) [kg] 
 

wA11 = gm * mA11w                # Gross weight of apollo 11 unar ascent vechicle [N] 

T11 = 16000                      # Thrust of apollo 11 lunar ascent vechicle [N] 
 

TtoW11 = T11 / wA11              # Thrust to weight ratio of apollo 11 lunar ascent vechicle 

minTtoW = 1.5                    # Assumed minimum Thrust to weight ratio of hypothetical lunar ascent vechicle w/ same thrust as apollo 11 one 
 

MassDiff = (T11 - wA11) / gm     # Mass needed to make Thrust to weight ratio of apollo 11 lunar ascent vechicle 0 [kg] 

# ^ mass of hypothetical lunar ascent vechicle cannot surpass this amount 

# since we are assume it will have the same thrust as the apollo 11 one 
 

massAvailable = ((T11/minTtoW)-wA11) / gm      # Available mass of lunar ascent vechicle to bring to lunar orbit [kg] 

# ^ mass maintains a 1.5 Thrust to weight ratio 


if (MassDiff - massAvailable) <= 0:  

    massAvailable = 0                          # Checks to make sure hypothetical mass does not exceed mass difference  
    

propWeightA11 = mA11w - mA11d     # Weight of propellant of apollo 11 lunar ascent vechicle [kg] 

propWeightD11 = mD11w - mD11d     # Weight of propellant of apollo 11 lunar descent vechicle [kg] 
 

IceBoxCost2a = mD11w * p                                   # Extra cost for The Ice Box per sale [$] 

IceBoxCost2b = mA11w * p                                   # Extra ONE TIME cost for The Ice Box [$] 

CustomerCost2a = 0                                         # Extra cost for customer per sale [$] 

CustomerCost2b = 0                                         # Extra ONE TIME cost for customer [$] 

SaleOFProductMass2 = massAvailable                         # Amount of product that can be transported per sale [kg] 

ReusableRequirements2 = propWeightA11                      # Amount of fuel needed per sale [kg] 
 

ROI = IceBoxCost2a + IceBoxCost2b 

MassPerHour = SaleOFProductMass2 / (((SaleOFProductMass2+ReusableRequirements2)/2000)*24) 

PowerPerSale = 0 

ParetoVector1 = [ROI, MassPerHour, PowerPerSale] 
 

print(ParetoVector1) 


#_________________________________________________________________________________________________________________ 

#_________________________________________________________________________________________________________________ 

## Falcon Heavy takes Apollo 11/Blue Origin type lander to Moon (Assumes only ascent vehicle is reusable) 
 

massAvailableF = 20000                     # FalconH mass to lunar surface [kg] 

priceFalcon = 100000000                    # Price to launch falcon heavy [$] 

p = priceFalcon / massAvailableF           # Price to launch one kg to the moon [$] (BASED on Falcon Heavy) 
 

# Apollo 11/Blue Origin type lander 

m11w = 32500 / lb2kg             # Wet mass of apollo 11 lunar lander (w/ propellant and crew) [kg] 

m11d = 9000 / lb2kg              # Dry mass of apollo 11 lunar lander (w/o propellant and crew) [kg] 
 

mA11w = 4700                     # Wet mass of apollo 11 lunar ascent vechicle (w/ propellant and crew) [kg] 

mA11d = 2150                     # Dry mass of apollo 11 lunar ascent vechicle (w/o propellant and crew) [kg] 
 

mD11w = m11w - mA11w             # Wet mass of apollo 11 lunar descent vechicle (w/ propellant and crew) [kg] 

mD11d = m11d - mA11d             # Dry mass of apollo 11 lunar descent vechicle (w/o propellant and crew) [kg] 
 

wA11 = gm * mA11w                # Gross weight of apollo 11 unar ascent vechicle [N] 

T11 = 16000                      # Thrust of apollo 11 lunar ascent vechicle [N] 
 

TtoW11 = T11 / wA11              # Thrust to weight ratio of apollo 11 lunar ascent vechicle 

minTtoW = 1.5                    # Assumed minimum Thrust to weight ratio of hypothetical lunar ascent vechicle w/ same thrust as apollo 11 one 
 

MassDiff = (T11 - wA11) / gm     # Mass needed to make Thrust to weight ratio of apollo 11 lunar ascent vechicle 0 [kg] 

# ^ mass of hypothetical lunar ascent vechicle cannot surpass this amount 

# since we are assume it will have the same thrust as the apollo 11 one  
 

massAvailable = ((T11/minTtoW)-wA11) / gm      # Available mass of lunar ascent vechicle to bring to lunar orbit [kg] 

# ^ mass maintains a 1.5 Thrust to weight ratio 


if (MassDiff - massAvailable) <= 0:  

    massAvailable = 0                          # Checks to make sure hypothetical mass does not exceed mass difference       

propWeightA11 = mA11w - mA11d     # Weight of propellant of apollo 11 lunar ascent vechicle [kg] 

propWeightD11 = mD11w - mD11d     # Weight of propellant of apollo 11 lunar descent vechicle [kg] 
 

IceBoxCost2a = mD11w * p                                   # Extra cost for The Ice Box per sale [$] 

IceBoxCost2b = mA11w * p                                   # Extra ONE TIME cost for The Ice Box [$] 

CustomerCost2a = 0                                         # Extra cost for customer per sale [$] 

CustomerCost2b = 0                                         # Extra ONE TIME cost for customer [$] 

SaleOFProductMass2 = massAvailable                         # Amount of product that can be transported per sale [kg] 

ReusableRequirements2 = propWeightA11                      # Amount of fuel needed per sale [kg] 


ROI = IceBoxCost2a + IceBoxCost2b 

MassPerHour = SaleOFProductMass2 / (((SaleOFProductMass2+ReusableRequirements2)/2000)*24) 

PowerPerSale = 0 

ParetoVector2 = [ROI, MassPerHour, PowerPerSale] 


print(ParetoVector2) 

#_________________________________________________________________________________________________________________ 

#_________________________________________________________________________________________________________________ 

## Starship takes Apollo 11/Blue Origin type lander to Moon (Assumes only ascent vehicle is reusable) 

massAvailableS = 100000                    # Starship mass to lunar surface [kg] 

priceStarship = 2000000                    # Price to launch falcon heavy [$] 

p = priceStarship / massAvailableS         # Price to launch one kg to the moon [$] (BASED on Starship)  

# Apollo 11/Blue Origin type lander 

m11w = 32500 / lb2kg             # Wet mass of apollo 11 lunar lander (w/ propellant and crew) [kg] 

m11d = 9000 / lb2kg              # Dry mass of apollo 11 lunar lander (w/o propellant and crew) [kg] 
 

mA11w = 4700                     # Wet mass of apollo 11 lunar ascent vechicle (w/ propellant and crew) [kg] 

mA11d = 2150                     # Dry mass of apollo 11 lunar ascent vechicle (w/o propellant and crew) [kg] 
 

mD11w = m11w - mA11w             # Wet mass of apollo 11 lunar descent vechicle (w/ propellant and crew) [kg] 

mD11d = m11d - mA11d             # Dry mass of apollo 11 lunar descent vechicle (w/o propellant and crew) [kg] 
 

wA11 = gm * mA11w                # Gross weight of apollo 11 unar ascent vechicle [N] 

T11 = 16000                      # Thrust of apollo 11 lunar ascent vechicle [N] 


TtoW11 = T11 / wA11              # Thrust to weight ratio of apollo 11 lunar ascent vechicle 

minTtoW = 1.5                    # Assumed minimum Thrust to weight ratio of hypothetical lunar ascent vechicle w/ same thrust as apollo 11 one 


MassDiff = (T11 - wA11) / gm     # Mass needed to make Thrust to weight ratio of apollo 11 lunar ascent vechicle 0 [kg] 

# ^ mass of hypothetical lunar ascent vechicle cannot surpass this amount 

# since we are assume it will have the same thrust as the apollo 11 one  
 
massAvailable = ((T11/minTtoW)-wA11) / gm      # Available mass of lunar ascent vechicle to bring to lunar orbit [kg] 

# ^ mass maintains a 1.5 Thrust to weight ratio 


if (MassDiff - massAvailable) <= 0:  

    massAvailable = 0                          # Checks to make sure hypothetical mass does not exceed mass difference  

     

propWeightA11 = mA11w - mA11d     # Weight of propellant of apollo 11 lunar ascent vechicle [kg] 

propWeightD11 = mD11w - mD11d     # Weight of propellant of apollo 11 lunar descent vechicle [kg] 


IceBoxCost2a = mD11w * p                                   # Extra cost for The Ice Box per sale [$] 

IceBoxCost2b = mA11w * p                                   # Extra ONE TIME cost for The Ice Box [$] 

CustomerCost2a = 0                                         # Extra cost for customer per sale [$] 

CustomerCost2b = 0                                         # Extra ONE TIME cost for customer [$] 

SaleOFProductMass2 = massAvailable                         # Amount of product that can be transported per sale [kg] 

ReusableRequirements2 = propWeightA11                      # Amount of fuel needed per sale [kg] 
 

ROI = IceBoxCost2a + IceBoxCost2b 

MassPerHour = SaleOFProductMass2 / (((SaleOFProductMass2+ReusableRequirements2)/2000)*24) 

PowerPerSale = 0 

ParetoVector3 = [ROI, MassPerHour, PowerPerSale] 

print(ParetoVector3)  

#_________________________________________________________________________________________________________________ 

#_________________________________________________________________________________________________________________ 

## Starship takes Starship type lander to Moon (Assumes vehicle is fully reusable) 

massAvailableS = 100000                    # Starship mass to lunar surface [kg] 

priceStarship = 2000000                    # Price to launch falcon heavy [$] 

p = priceStarship / massAvailableS         # Price to launch one kg to the moon [$] (BASED on Starship) 
 

# Starship type lander 

m11w = 2000000 / p                     # The variable name means nothing here this value is just meant to give cost of the  

                                       # ^starship vehicle for the print statement 
                                         

massAvailable = 100000                 # Available mass of lunar ascent/descent vechicle to bring to lunar orbit [kg] 

propWeightA11 = (1200*1000) / (g/gm)   # Weight of propellant of Starship vechicle needed per lunar launch [kg] 

propWeightD11 = 0                      # There is no descent vehicle so set propellant weight to 0 [kg] 
 

IceBoxCost2a = 0                                           # Extra cost for The Ice Box per sale [$] 

IceBoxCost2b = m11w * p                                    # Extra ONE TIME cost for The Ice Box [$] 

CustomerCost2a = 0                                         # Extra cost for customer per sale [$] 

CustomerCost2b = 0                                         # Extra ONE TIME cost for customer [$] 

SaleOFProductMass2 = massAvailable                         # Amount of product that can be transported per sale [kg] 

ReusableRequirements2 = propWeightA11 + propWeightD11      # Amount of fuel needed per sale [kg] 
 

ROI = IceBoxCost2a + IceBoxCost2b 

MassPerHour = SaleOFProductMass2 / (((SaleOFProductMass2+ReusableRequirements2)/2000)*24) 

PowerPerSale = 0 

ParetoVector4 = [ROI, MassPerHour, PowerPerSale] 


print(ParetoVector4) 

#_________________________________________________________________________________________________________________ 

#_________________________________________________________________________________________________________________ 

## Starship takes Spin Launch to Moon 

massAvailableS = 100000                    # Starship mass to lunar surface [kg] 

priceStarship = 2000000                    # Price to launch falcon heavy [$] 

p = priceStarship / massAvailableS         # Price to launch one kg to the moon [$] (BASED on Starship) 


spinw = 100000                                         # Weight of spin launch without vacuum infrastructure [kg] 

EarthMassLaunch = 200                                  # Mass of object launched on Earth [kg] 

MoonMassLaunch = (g/gm) * EarthMassLaunch              # Mass of object launched on Moon [kg] 

WeightOFequipment = 200                                # Weight of equipment needed to control product after launch [kg] 

V = 5000 / mph2ms                                      # Velocity of projectile at launch [m/s] 

Energy = (1/2) * MoonMassLaunch * V**2                  # Energy needed per launch [J]  

# Using Spin Launch 

IceBoxCost2a = 0                                           # Extra cost for The Ice Box per sale [$] 

IceBoxCost2b = spinw * p                                   # Extra ONE TIME cost for The Ice Box [$] 

CustomerCost2a = 0                                         # Extra cost for customer per sale [$] 

CustomerCost2b = 0                                         # Extra ONE TIME cost for customer [$] 

SaleOFProductMass2 = MoonMassLaunch - WeightOFequipment    # Amount of product that can be transported per sale [kg] 

ReusableRequirements2 = Energy                              # Amount of energy needed per sale [J] 


ROI = IceBoxCost2a + IceBoxCost2b 

MassPerHour = SaleOFProductMass2 / (((SaleOFProductMass2)/2000)*24) 

PowerPerSale = ReusableRequirements2 * (3/4)

ParetoVector5 = [ROI, MassPerHour, PowerPerSale] 


print(ParetoVector5)


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


def Power_Trades(Dist_Pwr, Dist_energy, Processing_Pwr, Processing_energy, a, b, c, d, e, f):

    #Useful variables
    cost_of_kg = 20 #Highly dependent on LV (this is based on starship values)
    sp_cost = 725/0.0032 #cost per square meter of solar panel

    Required_P = 28050 + Dist_Pwr + Processing_Pwr #(how much power the architecture needs, found from research estimates)
    Deg_mod = 1/(0.7**1.5) #Accounting for solar panel degredation (30% over ten years for 15 years)
    Power_systems = numpy.array([[0, Required_P*Deg_mod],[Required_P*0.3*Deg_mod/0.94,Required_P*0.7*Deg_mod],[Required_P*0.7*Deg_mod/0.94,Required_P*0.3*Deg_mod],[Required_P*Deg_mod/0.94,0]]) #Vector of hybrid architectures [Ground solar collection system, Orbital solar collection system]


    def orbit(E_req, x0, storage, D_R):
            D_T = 5 # Transmit diameter (m)
            wavelength = 0.00029979
            storage_cost, storage_mass, discharge_time = storage

            eff, rectenna_cost = orbit_power_output(x0, D_T, D_R, wavelength, discharge_time)
            
            TU_to_s = 382981
            P_req = E_req / (x0[-1]-discharge_time*3600/TU_to_s)
            P_out = P_req / eff

            num_satellites = cal_satellites(P_out)

            A_solar = 200

            kgpa = 7.8125
            m_solar = A_solar * kgpa

            cost_of_kg = 20
            baseline_sat_cost = 50_000_000
            baseline_sat_mass = 1000 # 5100 going into transfer orbit

            sp_cost = 725/0.0032 #cost per square meter of solar panel

            total_mass = baseline_sat_mass + m_solar + storage_cost
            sat_cost = baseline_sat_cost + sp_cost * A_solar + total_mass * cost_of_kg + storage_cost

            A_R = (D_R/2)**2 * math.pi # Reciever aperature (m^2)
            rectenna_mass = A_R * sp_cost
            rectenna_cost += rectenna_mass * cost_of_kg

            power = P_out
            fuel = 0
            Cost = sat_cost * num_satellites + rectenna_cost
            lifetime = 10

            return power, fuel, Cost, lifetime


    def beaming(D_T, D_R, d, wavelength):
        A_T = (D_T/2)**2 * math.pi # Transmit aperature (m^2) 
        A_R = (D_R/2)**2 * math.pi # Reciever aperature (m^2)

        rectenna_cost_per_m2 = 460
        rectenna_cost = rectenna_cost_per_m2 *A_R
        
        tau = (A_R*A_T)**0.5/(wavelength*d) # Efficiency parameter for design parameters
        eff = (1-math.e**(-tau**2) )# Design Efficiency

        return eff, rectenna_cost

    def sat_batt(E_req):
        # Constants
        energy_density = 70 # Wh/kg
        const_current = 500 # A discharging
        pulse_current = 2000 # A discharging
        V = 3.6 # V
        discharge_percent = 0.3
        cost_per_kWh = 1000 # $
        energy = energy_density * 0.275 * discharge_percent  # W h

        power = pulse_current*V  # W 
        num = E_req / energy
        mass = 0.275 * num
        
        cost = cost_per_kWh * energy *1000
        discharge_time = energy / power

        return cost, mass, discharge_time

    def sat_supercap(E_req):
        # Constants
        energy_density = 5.69 # Wh/kg
        power_density = 7916.66 # W/kg
        V = 2.85 # V
        cost_per_kWh = 11_792
        discharge_percent = 0.98

        energy = energy_density*0.072 * discharge_percent
        power = power_density*0.072
        num = E_req /energy
        mass = 0.072 * num
        
        discharge_time = energy / power
        cost = cost_per_kWh * energy *1000

        return cost, mass, discharge_time

    def orbit_power_output(x0, D_T, D_R, wavelength, discharge_time):
        sys = System(P1=Earth, P2=Moon)
        x1, P1 = np.split(x0,[-1])

        ode_kw = dict(rel_tol=1E-13, abs_tol=1E-13) # integration parameters
        sol1 = cr3bp.propagate(sys, x1, np.array([0, 10*P1]), **ode_kw)

        min_dist = 100000000000000000
        LU_to_KM = 389703
        moon_radius = 1737.1
        min_idx = -1
        for i in range(len(sol1.y[0])):

            dist = np.sqrt((LU_to_KM*sol1.y[0,i])**2 + (LU_to_KM*sol1.y[1,i])**2 + (sol1.y[2,i]*LU_to_KM-moon_radius)**2)
            if dist < min_dist:
                min_idx = i
                min_dist = dist


        end_discharge_idx = min_idx
        TU_to_s = 382981
        while sol1.t[end_discharge_idx] < (sol1.t[min_idx] + discharge_time*3600/TU_to_s):
            end_discharge_idx += 1

        end_dist = (LU_to_KM*sol1.y[0,end_discharge_idx])**2 + (LU_to_KM*sol1.y[1,end_discharge_idx])**2 + (sol1.y[2,end_discharge_idx]*LU_to_KM-moon_radius)**2

        eff, rectenna_cost = beaming(D_T, D_R, min_dist, wavelength)

        return eff, rectenna_cost

    def cal_satellites(Picebox):
        # determine Psa, how much power solar array needs to provide during daylight periods
        Pe = 0 # [W] power required by sc during eclipse
        Te = 0 # [s] period in eclipse
        Xe = 0.65 # [] efficiency of pasth from Solar array through batteries to loads
        Pd = 1000 # [W] power required by sc during daylight periods
        Td = 635040 # [s] period in daylight, assuming always in daylight
        #T = np.array([635040, 159408, 19926, 2551392]) # [sec]
        Xd = 0.80 # [] efficiency of pasth from Solar array to loads directly

        Psa = (Pe*Te/Xe + Pd*Td/Xd) / Td

        # determine Pbol, power reqd at beginning of life
        n_solar = 0.4 # solar efficiency
        I = 1365 # [W/m^2] solar constant
        P0 = I * n_solar # [W/m^2], ideal solar cell output performance
        Id = 0.77 # [], nominal solar array degredation
        theta = 0 #deg, assumes solar incidence is normal throughout orbits

        Pbol = P0*Id*1 #np.cos(theta) #[W/m^2]

        # determine Peol, power reqd at end of life
        Ld = 0.75 #[%] degredation assuming 15 year life degredation, usually due to radiation
        Peol = Pbol*Ld # [W/m^2]

        # solar array area needed to power spacecraft
        Asa = Psa/Peol # [m^2]
        # NOTE: main variables are time on orbit (in sun), and power req'd by S/C to function

        # Now make assumption that we can double our solar arrays to 
        # VARY this in a trade up to 75 m^2
        # Max is 400m^2
        Atot = 200 #[m^2], total Area on spacecraft
        Apwr = Atot - Asa #[m^2] how much area of solar panels we can use for power transmission

        # How much power we can supply by each satellite
        Seff = 0.15 #typical solar panel efficiency
        Beff = 1 # power losses to beam power
        SI = 1366.1 # [W/m^2] solar irradiance near L1
        Pt = Apwr*Seff*Beff*SI # [W] how much power we can beam, assuming 50% losses

        # Determine how many satellites we need
        # Picebox = 1*(10**5) #[W]
        numSats = Picebox/Pt
        
        return numSats
    
    ##POWER TRADES

    # Input array
    orbit_arr = np.array([
    # Halo
    [1.03371467823365E+00,	0.0000000000000000E+00,	1.89133683287009E-01,	-3.95740371598041E-14,	-1.27360767906067E-01,	8.55682270827113E-13,	1.6663353286308E+00],
    # Butterfly
    [9.1202100067489833E-1,	0.0000000000000000E+0,	1.4950993789808714E-1,	1.5379452260107897E-13,	-2.7261732734822092E-2,	1.1488748710526696E-12,	3.3282518875369433E+0],
    # Distant Prograde
    [1.0236206299055113E+0,	0.0000000000000000E+0,	0.0000000000000000E+0,	5.5491759614030665E-14,	5.5040444235154051E-1,	0.0000000000000000E+0,	4.1636825485004247E-1],
    # Lyapunov
    [9.9690781458801814E-1,	0.0000000000000000E+0,	0.0000000000000000E+0,	-6.9583985002249150E-14,	1.6447667663718337E+0,	0.0000000000000000E+0,	6.6624743663890538E+0]
    ])

    #Ground solar power system
    Ground_Solar = Power_systems[b,0]
    #print(Ground_Solar)

    n_solar = 0.4
    I = 1365
    A_ground = Ground_Solar/(I*n_solar) #Area of ground solar panels
    #print(A_ground)
    kgpa = 7.8125
    m_ground = A_ground * kgpa
    Ground_Solar_Lifetime = 15

    #Ground Batteries

    Energy_Margin = Required_P*60*60*24*3
    Stockpile_Pwr = Dist_Pwr + Processing_Pwr
    Storage_req = ((Dist_energy + Processing_energy + Energy_Margin)/3600) #Wh required for max storage
    #print(Storage_req)

    if a == 0:
        #1: Deep impact
        SE_DI = 250 #Specific energy of deep impact battery Wh/kg
        mass_DI = Storage_req/SE_DI
        #print(mass_DI)

        Battery_Cost = mass_DI*cost_of_kg
        Battery_Lifetime = 15


    elif a == 1:
        #2: Flywheels
        SE_FW = 2175480/3600
        mass_FW = Storage_req/SE_FW

        Battery_Cost = mass_FW * cost_of_kg
        Battery_Lifetime = 20


    #Orbital solar power system
    Orbital_Solar_Req = Power_systems[b,1]
    D_R = [20,25,30] # m





    if e == 0:
        Orbit_Pwr, Orbit_Fuel, Orbit_Cost, Orbital_Lifetime = orbit(Orbital_Solar_Req,orbit_arr[d],sat_batt(Required_P),D_R[f])
    elif e == 1:
        Orbit_Pwr, Orbit_Fuel, Orbit_Cost, Orbital_Lifetime = orbit(Orbital_Solar_Req,orbit_arr[d],sat_supercap(Required_P),D_R[f])



    #RTG's
    RTG_power = numpy.array([0,110])


    GS_cost = A_ground*sp_cost + m_ground*cost_of_kg #Cost of ground solar sytsem
    #print(GS_cost)
    RTG_costs = 100000000*numpy.array([0,1,2,3])#Cost of RTG's
    RTG_Cost = RTG_costs[c]



    
    #print("\n\n",e,"\n\n")

    Power = Required_P + RTG_power[c]
    Fuel = 0
    Cost = GS_cost + Battery_Cost + RTG_Cost + Orbit_Cost
    Efficiency = Power
    Lifetime = min(Battery_Lifetime,Ground_Solar_Lifetime,Orbital_Lifetime)

    return (Power,Fuel,Cost,Efficiency,Lifetime)

# return power, fuel, cost, efficiency, lifetime


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

Dist_Pwr = ParetoVector5[2]/(12*3600)
Dist_energy = ParetoVector5[2]
a = [0,1]
b = [0,1,2,3]
c = [0,1]
d = [0,1,2,3]
e = [0,1]
f = [0,1,2]

#Power,Fuel,Cost,Efficiency,Lifetime

counter = 0

loops = len(a)*len(b)*len(c)*len(d)*len(e)*len(f)

#Doug = np.zeros((loops,1))
#John = Doug
#Bill = John
#Fred = Bill
#George = Fred

ProcessingPower = [ProcBuckPower[0],ProcBuckPower[1],ProcConvPower[0],ProcConvPower[1],ProcDrumPower[0],ProcDrumPower[0]]

ProcessingEnergy = ProcessingPower * (12*3600)
header = ['Power','Fuel','Cost','Efficiency','Lifetime','Power Code']
## Loop through power model and comine all options

with open('power.csv', 'w') as file:
    writer = csv.writer(file)
    writer.writerow(header)

    for q in range(0,1): #Power
        Processing_Pwr = ProcessingPower[q]
        Processing_energy = ProcessingEnergy[q]
        for i in range(0,1): #len(a)): #Battery
            a1 = a[i]
            #print(i)
            for j in range(0,4): # Hybrid Architecture ONLY ONE THAT CHANGES INFRASTRUCTURE COST
                b1 = b[j]
                for k in range(0,1): # RTG's
                    c1 = c[k]
                    for z in range(0,1): #Orbits
                        d1 = d[z]
                        for y in range(0,1): # Orbit Battery
                            e1 = e[y]
                            for x in range(0,1): # Reciever Diameter
                                f1 = f[x]
                                print("\nThis is loop number:",counter,"/",loops,"\n Values for loop PowerModel:",Processing_Pwr," a:",a1," b:",b1," c:",c1," d:",d1," e:",e1," f:",f1,"\n\n")
                                [Doug,John,Bill,Fred,George] = Power_Trades(Dist_Pwr, Dist_energy, Processing_Pwr, Processing_energy, a1, b1, c1, d1, e1, f1)
                                #print(Doug,John,Bill,Fred,George)
                                PowerCode = str((q,a1,b1,c1,d1,e1,f1))
                                writer.writerow([Doug,John,Bill,Fred,George,PowerCode])
                                counter = counter + 1
                                #print("\n",i," ",j," ",k," ",z," ",y," ",x,"\n")

                            if counter == 384:
                                counter = 0
                                writer.writerow(["-","-","-","-","-","-"])



                            x = x+1
                        y = y+1
                    z = z+1
                k = k+1
            j = j+1
        i = i+1
    q = q + 1

file.close()

#print("\n","\n",Doug,"\n",John,"\n",Bill,"\n",Fred,"\n",George,"\n")


##ROI TRADES -----------------------------------------------------------------------------

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
OS_Eff=100; #kW
OS_Cost=110000000; #USD

OS_USD_per_kJ=OS_Cost*(1/OS_Eff)*(1/OS_LS)*(1/31536000); #$/kJ


#Ground Solar Array


GS_LS=12;  #Years
GS_Eff=100; #kW
GS_Cost=90000000; #USD

GS_USD_per_kJ=GS_Cost*(1/GS_Eff)*(1/GS_LS)*(1/31536000); #$/kJ





#70/30 Orbital Ground Mix

OG70_LS=14;  #Years
OG70_Eff=100; #kW
OG70_Cost=105000000; #USD

OG70_USD_per_kJ=OG70_Cost*(1/OG70_Eff)*(1/OG70_LS)*(1/31536000); #$/kJ


#30/70 Orbital Ground Mix

OG30_LS=13;  #Years
OG30_Eff=100; #kW
OG30_Cost=95000000; #USD

OG30_USD_per_kJ=OG30_Cost*(1/OG30_Eff)*(1/OG30_LS)*(1/31536000); #$/kJ


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



Ground_Eff=[ProcMethBucket[0],ProcHydrBucket[0],ProcMethDrum[0],ProcHydrDrum[0],ProcMethConv[0],ProcHydrConv[0]] 
Ground_LS=[BL_LS,BL_LS,BC_LS,BC_LS,DL_LS,DL_LS]
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


Revenue = np.zeros((6,7))


#print(float(Ground_Eff[0])*float(Ground_LS[0])*(870)*float(Product_Abd[3])*float(Product_Cost[3]))

count = 1

for i in range(0,6): #Ground
    for j in range(0,7): #Products
        #print(count)
        #count = count + 1
        #print(float(Ground_Eff[i])*float(Ground_LS[i])*(870)*float(Product_Abd[j])*float(Product_Cost[j]))
        Revenue[i,j] = float(Ground_Eff[i])*float(Ground_LS[i])*(870)*float(Product_Abd[j])*float(Product_Cost[j]);

#print("Revenue")
#print(Revenue)

Power_Eff=[GS_Eff,OG30_Eff,OG70_Eff,OS_Eff];
Power_LS=[GS_LS,OG30_LS,OG70_LS,OS_LS];

Power_Cost=[GS_USD_per_kJ,OG30_USD_per_kJ,OG70_USD_per_kJ,OS_USD_per_kJ,];
Ground_Cost=[BL_USD_per_kg,BL_USD_per_kg,BC_USD_per_kg,BC_USD_per_kg,DL_USD_per_kg,DL_USD_per_kg];
Processing_Cost=[CR_USD_per_kg,HR_USD_per_kg];
Transportation_Cost=[T1_USD_per_kg,T2_USD_per_kg,T3_USD_per_kg,T4_USD_per_kg,T5_USD_per_kg];

#Infastructure_Cost= np.zeros((6,30))
Infastructure_Cost= np.zeros((len(Ground_Eff),len(Power_Cost)*len(Ground_Cost)*len(Transportation_Cost)))

#print(Power_Eff)
#print(Power_LS)
#print(Power_Cost)
#print(Ground_Cost)
#print(Processing_Cost)
#print(Transportation_Cost)

ROICode = []

overall_count = 0

k= 0;
for i in range(0,6): #Ground
    counter=0
    #print("i:",i)
    for j in range(0,4): #Power
        #print("j:",j)
        for h in range(0,5): #Transportation
            
            if k == 2: # if even PROCESSING METHOD
                k = 0

            #print("h:",h)
            #print("k:", k)
            Infastructure_Cost[i,counter] = float(Ground_LS[i])*float(Ground_Eff[i])*870*(float(Ground_Cost[i])+float(Processing_Cost[k])+float(Transportation_Cost[h]))+float(Power_LS[j])*float(Power_Eff[j])*31536000*float(Power_Cost[j])
            ROICode.append(str((i,j,h,k)))
            counter += 1 
            overall_count += 1
 
    k += 1


print(ROICode)



#print("Infastructure_Cost")
#print(Infastructure_Cost)

#Only 3 Ground Options

#Each Individual Ground Option has 7 Revenue Options, and 12 infastructure
#cost options. Need to Compare each Revenue Option to 12 Infastrcture Cost
#Options 

Revenue_vs_Cost = np.zeros((6*7*12,2)) #[[0 for col in range(3*7*12)] for row in range(2)]
counter = 0
for i in range(0,6): #Mining Option
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
#print("Revenue v Cost")       
#print(Revenue_vs_Cost)

Profit= Revenue_vs_Cost[:,0] - Revenue_vs_Cost[:,1]
print("Profit")
print(Profit)


        #RoiScore = Profit
        #PowerScore = a b c d e f
        #FuelScore = of the -6 (Bucket, Drum, Conv) (Methane, Hydrogen)

        #[Fuel, Power, ROI,ROIScore,PowerScore,FuelScore]



#---------------------------------------------------------------------------------------------
#power initial costs
#GS_cost = A_ground*sp_cost + m_ground*cost_of_kg #Cost of ground solar sytsem
#print(GS_cost)
#OS_cost = m_orbital*cost_of_kg + num_of_sat*cost_of_sat#Cost of orbital solar system
RTG_cost = 100000000*numpy.array([0,1,2,3])#Cost of RTG's


##Final Pareto
count = 1

#for a = 1:
#    for b = 1:
#        for c = 1:
#            for d = 1:
#                for e = 1:
#                    for f = 1:
#                        Power(count) = 
#                        Fuel(count) = 
#                        ROI(count) = 
#                        count = count +1
#

#Plot Pareto from the three vectors above
# Fixing random state for reproducibility
np.random.seed(19680801)


def randrange(n, vmin, vmax):
    """
    Helper function to make an array of random numbers having shape (n, )
    with each number distributed Uniform(vmin, vmax).
    """
    return (vmax - vmin)*np.random.rand(n) + vmin

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

n = 100

# For each set of style and range settings, plot n random points in the box
# defined by x in [23, 32], y in [0, 100], z in [zlow, zhigh].
for m, zlow, zhigh in [('o', -50, -25), ('^', -30, -5)]:
    xs = randrange(n, 23, 32)
    ys = randrange(n, 0, 100)
    zs = randrange(n, zlow, zhigh)
    ax.scatter(xs, ys, zs, marker=m)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()