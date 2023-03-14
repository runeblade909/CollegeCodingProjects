#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov  8 21:04:14 2022

@author: williamgettinger
"""

import numpy as np

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
P0 = 190 # [W/m^2], lower bound of ideal solar cell output performance
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
Picebox = 1*(10**6) #[W]
numSats = Picebox/Pt

# 2 batteries, beaming and self-sustain
# determine size of beaming battery based off of Pt

# Battery Parameters - using Deep Impact LI-SOCl2    
V = 24 # [volts] min value for now...
C = 312 # [Ah]
m = 36.6 # [kg]
e = 250 # [Wh/kg]
dl = 6 # [yr] ... incorperate into Pareto later (i.e: ROI)





