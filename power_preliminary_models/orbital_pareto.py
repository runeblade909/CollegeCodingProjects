from csltk import cr3bp
from csltk.utilities import Earth, Moon, System
from csltk.correctors import SingleShooter

import numpy as np
import matplotlib.pyplot as plt
import math

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

    kgpa = 321.16
    m_solar = A_solar * kgpa

    cost_of_kg = 20
    baseline_sat_cost = 50_000_000
    baseline_sat_mass = 1000 # 5100 going into transfer orbit

    sp_cost = 500/0.0032 #cost per square meter of solar panel

    total_mass = baseline_sat_mass + m_solar + storage_cost
    sat_cost = baseline_sat_cost + sp_cost * A_solar + total_mass * cost_of_kg + storage_cost

    A_R = (D_R/2)**2 * math.pi # Reciever aperature (m^2)
    rectenna_mass = A_R * sp_cost
    rectenna_cost += rectenna_mass * cost_of_kg

    power = P_out
    fuel = 0
    ROI = sat_cost * num_satellites + rectenna_cost
    lifetime = 10

    return power, fuel, ROI, lifetime

print(orbit(214000, orbit_arr[0], sat_batt(214000),20))
print(orbit(214000, orbit_arr[0], sat_supercap(214000),20))

# Input arrays
# D_R = [20,25,30] # m
# storage = [sat_batt(E_req),sat_supercap(E_req)]