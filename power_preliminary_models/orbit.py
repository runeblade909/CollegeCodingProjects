from csltk import cr3bp
from csltk.utilities import Earth, Moon, System
from csltk.correctors import SingleShooter

import numpy as np
import matplotlib.pyplot as plt
import math

######## TODO
##### Account for battry and solar array degredation

def solar_array(area):
    return

def sat_batt(num):
    # Constants
    energy_density = 70 # Wh/kg
    const_current = 500 # A discharging
    pulse_current = 2000 # A discharging
    V = 3.6 # V
    discharge_percent = 0.5

    mass = 0.275 * num


    power = pulse_current*V/2  # W 
    energy = energy_density * mass * discharge_percent  # W h

    discharge_time = energy / power

    return energy, power, mass, discharge_time


def sat_supercap(num):
    # Constants
    energy_density = 5.69 # Wh/kg
    power_density = 7916.66 # W/kg
    V = 2.85 # V
    cost_per_kWh = 11_792
    discharge_percent = 0.98


    mass = 0.072 * num
    power = power_density*mass
    energy = energy_density*mass * discharge_percent

    discharge_time = energy / power

    return power, energy, mass, discharge_time


E_out_batt, P_out_batt, batt_mass, discharge_time = sat_batt(100)

# define Earth-Moon system and orbit parameters
sys = System(P1=Earth, P2=Moon)
# ax = sys.plot_system()

# Halo
x0_guess = np.array([1.03371467823365E+00,	0.0000000000000000E+00,	1.89133683287009E-01,	-3.95740371598041E-14,	-1.27360767906067E-01,	8.55682270827113E-13,	1.6663353286308E+00])

# Butterfly
# x0_guess = np.array([9.1202100067489833E-1,	0.0000000000000000E+0,	1.4950993789808714E-1,	1.5379452260107897E-13,	-2.7261732734822092E-2,	1.1488748710526696E-12,	3.3282518875369433E+0])

# Distant Prograde
x0_guess = np.array([1.0236206299055113E+0,	0.0000000000000000E+0,	0.0000000000000000E+0,	5.5491759614030665E-14,	5.5040444235154051E-1,	0.0000000000000000E+0,	4.1636825485004247E-1])

# Lyapunov
# x0_guess = np.array([9.9690781458801814E-1,	0.0000000000000000E+0,	0.0000000000000000E+0,	-6.9583985002249150E-14,	1.6447667663718337E+0,	0.0000000000000000E+0,	6.6624743663890538E+0])

x1, P1 = np.split(x0_guess,[-1])

# x1,P1,_ = SingleShooter.fixed_time(x0_guess, sys, tol=1e-10) # period is fixed
# x1,P1,_ = SingleShooter.variable_time(x0_guess, sys, tol=1e-10) # period is free to vary

# Orbit solution, simulation of 1 period
# Units in LU and TU from API but converted down below
ode_kw = dict(rel_tol=1E-13, abs_tol=1E-13) # integration parameters
sol1 = cr3bp.propagate(sys, x1, np.array([0, 50*P1]), **ode_kw)


print('Orbit initial conidition =',x1)
print('Orbit period =',P1)

# ax.plot(sol1.y[0,:], sol1.y[1,:], sol1.y[2,:], label='Corrected')
# ax.legend()

# plt.show()

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
while sol1.t[end_discharge_idx] < (sol1.t[min_idx] + discharge_time*3600/TU_to_s/2):
    end_discharge_idx += 1

###### TODO    
###### Adjust efficiency for discharge time #####

print('Time of min distance =', sol1.t[min_idx],' Min distance =', min_dist, 'km')
discharge_idx_diff = end_discharge_idx-min_idx
start_idx = min_idx - discharge_idx_diff
print('Beaming start and stop times (days) =', sol1.t[min_idx-discharge_idx_diff]*TU_to_s/3600/24, sol1.t[end_discharge_idx]*TU_to_s/3600/24, sol1.t[-1]*TU_to_s/3600/24)

end_dist = (LU_to_KM*sol1.y[0,start_idx])**2 + (LU_to_KM*sol1.y[1,start_idx])**2 + (sol1.y[2,start_idx]*LU_to_KM-moon_radius)**2
# incidence = math.acos(end_dist/min_dist)

print('incidence =',min_dist,end_dist)

# Design Parameters 
D_T = 4.8768 # Transmit diameter (m)
A_T = (D_T/2)**2 * math.pi # Transmit aperature (m^2) 
D_R = 20 # Reciever diameter (m)
A_R = (D_R)**2  # Reciever aperature (m^2)
A_R = (D_R/2)**2 * math.pi

d = min_dist # transmit distance (m)
freq = 700 # Frequency (GHz)
wavelength = 0.00029979

tau = (A_R*A_T)**0.5/(wavelength*d) # Efficiency parameter for design parameters
eff = (1-math.e**(-tau**2) )# Design Efficiency


print('Beaming efficiency =',eff)


P_out = P_out_batt * eff
E_out = E_out_batt * eff


###############################
# Normalized for pareto
# Need battery and solar array worked out

P_out = P_out * 100
E_out = E_out * 100

print('Power output =',P_out, 'W')
print('Energy output =',E_out,'Wh')

baseline_sat_cost = 100_000_000
baseline_sat_mass = 1883 # 5100 going into transfer orbit
baseline_antenna_area = (11/2)**2 * math.pi # \ft^2
baseline_solar_array_area = 2 * (13.2**2) # m^2
# solar_array_cost_per_m2 = 
# solar_array_mass_per_m2 = 
# transmit_dish_cost_per_m2 =
# tansmit_dist_mass_per_m2 =

cost_of_kg = 20000

############ TODO
############ Account for cost of specific launch trajectories

total_mass = baseline_sat_mass + batt_mass
sat_cost = baseline_sat_cost + total_mass * cost_of_kg

print('Satellite cost including launch =',sat_cost)

power_solar = 400
n_solar = 0.4
I = 1365
A_solar = power_solar/(I*n_solar)

kgpa = 321.16
m_solar = A_solar * kgpa

print(power_solar)