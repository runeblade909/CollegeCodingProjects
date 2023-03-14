import math

# Input Parameters from past mission
# Current inputs are NASA's reference system parameters
# D_T1 = 1000 * 10 **-3# Transmit diameter (m)
# A_T1 = (D_T1/2)**2 * math.pi # Transmit aperature (m^2) 
# D_R1 = 10000 * 10 ** -3# Reciever diameter (m)
# A_R1 = (D_R1/2)**2 * math.pi # Reciever aperature (m^2)
# P_out1 = 8500 # Power output (MW)
# eff_1 = 0.588 # DC-DC efficiency
# d_1 = 35786 # transmit distance (km)
# freq1 = 2.45 # Frequency (GHz)


D_T1 = 1000 * 10 **-3# Transmit diameter (m)
A_T1 = 30*30 # Transmit aperature (m^2) 
D_R1 = 10000 * 10 ** -3# Reciever diameter (m)
A_R1 = 1000*1000 # Reciever aperature (m^2)
P_out1 = 10 # Power output (MW)
eff_1 = 0.588 # DC-DC efficiency
d_1 = 400_000 # transmit distance (km)
freq1 = 10 # Frequency (GHz)
wavelength_1 = 0.02997925 


# Design Parameters 
D_T2 = 4.8768 # Transmit diameter (m)
A_T2 = (D_T2/2)**2 * math.pi # Transmit aperature (m^2) 
D_R2 = 40 # Reciever diameter (m)
A_R2 = (D_R2)**2  # Reciever aperature (m^2)
P_out2 = 1 # Power output (MW)
d_2 = 4_250_000 - 3_475_000 # transmit distance (m)
freq2 = 700 # Frequency (GHz)
wavelength_2 = 0.00029979



tau_1 = (A_R1*A_T1)**0.5/(wavelength_1*d_1) # Efficiency Parameter for reference parameters
tau_2 = (A_R2*A_T2)**0.5/(wavelength_2*d_2) # Efficiency parameter for design parameters

eff_1 = (1-math.e**(-tau_1**2) )
eff_2 = (1-math.e**(-tau_2**2) )# Design Efficiency

print(tau_1,tau_2)
print(eff_1)
print(eff_2)


