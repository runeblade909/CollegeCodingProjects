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
import pandas as pd
import statistics
import math
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.cm as cm


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


	


