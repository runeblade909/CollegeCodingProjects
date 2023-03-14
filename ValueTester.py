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


"""
Processing Code
    Fabrizio Roberts

This code compares two regolith processing methods: Carbothermal and Hydrogen Reduction

%% Assumptions:
%   - starting from T= 260K
%   - ideal stoichiemtry
%   - ideal reaction conditions
%   - No slag or excess CO2 produced
%   - Perfect scaling
%   - Adiabatic container
%   - No surface film

"""
from math import log,exp,pow
import numpy as np

#Processing functions


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
    Reg = 1000; # Regolith (kg)
    scale = regRate/Reg; 
    m = Reg; # mass of regolith 
    
    Tday = 260; #lunar surface temp at equator during day
    Tnight = 140; #lunar surface temp at equator during night
    time = Reg/(regRate*3600); #time since start of process
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
    
    Ptot = Ptot * scale;
    #time = time * scale;
    
    Metals = [(Ti*molTi)*scale,(Si*molSi)*scale,(Mg*molMg)*scale,(Al*molAl)*scale,(Fe*molFe)*scale,(Ca*molCa)*scale]; #kg
    Power = [Ptot/time,0]; #J/s
    Time = [time/3600,0]; #hr
    Cost = [scale*(CH4req*pPkg*(molC+4*molH)),0]; #$
    Fuel = [scale*(H2tot*2*molH) + O2tot*2*molO, 0]; #kg
    Weight = [np.sum(Metals[0])+Fuel[0],0]; #kg
    
    """
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
    """

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
    
    Ptot = Ptot * scale;
    #time = time * scale;
    
    #Metals[1] = [(Ti*molTi)*scale,(Si*molSi)*scale,(Mg*molMg)*scale,(Al*molAl)*scale,(Fe*molFe)*scale,(Ca*molCa)*scale]; #kg
    Power[1] = Ptot/time; #J/s
    Time[1] = time/3600; #hr for total weight
    Cost[1] = scale*(H2req*(2*molH)*pPkg); #$
    Fuel[1] = scale*((H2tot*2*molH) + O2tot*2*molO); #kg
    Weight[1] = np.sum(Metals[1])+Fuel[1]; #kg
  

    """
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
    """
    # hr/kg over 1000 kg

    Tt1000 = np.divide(Time,Weight) * 1000
    Fuel = np.divide(Fuel,Time) #kg -> kg/hr

    return(Metals,Power,Cost,Fuel,Time,Tt1000,scale)


ProcBuckMetal,ProcBuckPower,ProcCost,ProcBuckFuel,ProcTime,ProcBuckTime1000,scale = Processing(16200)


print("\n",ProcTime,"\n",ProcBuckPower,"\n")