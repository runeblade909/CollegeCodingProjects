

from cmath import pi

import numpy as np

## This code is to evaluate

#Excavation Rate
#Power Consumption
#Simplicity

class ExcavationAnalysis:
	def __init__(self):
		#Regolith Properties
		self.Rho = 1.5 #g/cm^3
		self.Rho = self.Rho *100**3/1000 #Kg.m^3
		#Power Constraints
		self.MotorPower = 50 # Each motor takes 50W
		self.TorqueMotor = 75 # Higher torque application motor

		#Bucket loaded
		self.VolBucket = .6*.5*.5  # Volume of bucket m^3

		self.ExcMassBucket = self.Rho * self.VolBucket  #Bucket mass 

		self.ExcRateBucket = self.ExcMassBucket *.75 * 60#Kg/hr

		self.SimplicityBucket = 8

		self.PowerBucket = self.TorqueMotor*2 #2 motor linear actuators to raise and move bucket 

		#Drum Excavator
		self.r = .25 #m
		self.L = 0.5 #m
		self.VolDrum = pi*self.r**2 *self.L #m^3 drum volume

		self.ExcMassDrum = self.Rho * self.VolDrum  #kg

		self.ScoopSize = .1*.1*.2   #m^3 for 5 scoops across the length, 3 distinct rows (see NASA Rassor)

		self.ExcRateDrum = 12 * (self.ScoopSize*15)*self.Rho*60 # Kg/hr

		self.PowerDrum = self.TorqueMotor*4 #2 motors to move bucket and rotate bucket

		#Bucket Conveyor
		self.OneBucket = .5*.1*.1

		self.VolBin = .7*.3*.5 #Volume of storage m^3 bin for conveyor belt regolith

		self.MassBin = self.VolBin*self.Rho # Mass kg

		self.VolConv = self.OneBucket * 10

		self.ExcMassConv = self.Rho * self.VolConv/2

		self.ExcRateConv = self.ExcMassConv * 12 * 60 #Kg/hr

		self.PowerConv = self.MotorPower + 4*self.TorqueMotor # 1 linear actuator motor to lower to ground, 2 motors for rotating coveyor and positioning

	def excavator_productivity(self):
		#I used this resource for productivity
		#https://homesteady.com/13650496/how-to-calculate-excavator-productivity

		#Q = (60*q*z*n*kf)/kl [cubic ft/hr]
		#q = volume of bucket
		#z = number of buckets
		#n = RPM or speed of rotation of the rotor
		#kf = filling factor of bucket (0-1)
		#kl = soil-loosening factor ~1.2

		self.QBucket = (60*self.VolBucket*1*1.5*1*.7)/1.2
		self.QDrum = (60*self.ScoopSize*15*12*.85)/1.2
		self.QConv = (60*self.OneBucket*10*12*.8)/1.2

		self.RateBucket = self.QBucket * self.Rho
		self.RateDrum = self.QDrum * self.Rho
		self.RateConv = self.QConv * self.Rho

		## Final vector scores
		self.BucketVec = np.array([self.RateBucket,self.PowerBucket,0])
		self.DrumVec = np.array([self.RateDrum,self.PowerDrum,0])
		self.ConvVec = np.array([self.RateConv,self.PowerConv,0])

if __name__ == '__main__':

	ExcavationAnalysis = ExcavationAnalysis()
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


	print(ExcavationAnalysis.ScoopSize)
