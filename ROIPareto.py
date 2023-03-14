# Pareto Code
# By: Amir Tillis, Zak, Fabrizio
# Last Updated: 11/17/22
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

def ROI(MethBuckF,HydrBuckF,MethDrumF,HydrDrumF,MethConvF,HydrConvF,BuckMetal,DrumMetal,ConvMetal,scale,ProcCost,ProcTime):


    #Orbital Solar Array b=3

    OS_LS=10;  #Years
    OS_Eff=100; #W
    OS_Cost=110000000; #USD

    OS_USD_per_kJ=OS_Cost*(1/OS_Eff)*(1/OS_LS)*(1/31536000); #$/kJ


    #Ground Solar Array b = 0


    GS_LS=15;  #Years
    GS_Eff=100; #kW
    GS_Cost=50000000; #USD

    GS_USD_per_kJ=GS_Cost*(1/GS_Eff)*(1/GS_LS)*(1/31536000); #$/kJ


    # 30/70 Ground Solar b = 1
    OG30_LS = 11.5
    OG30_Eff = 100
    OG30_Cost = 127000000

    OG30_USD_per_kJ=OG30_Cost*(1/OG30_Eff)*(1/OG30_LS)*(1/31536000); #$/kJ

    #70/30 Orbital Ground Mix b = 2

    OG70_LS = 13.5
    OG70_Eff = 100
    OG70_Cost = 83000000

    OG70_USD_per_kJ=OG70_Cost*(1/OG70_Eff)*(1/OG70_LS)*(1/31536000); #$/kJ
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
    CR_Eff=np.multiply(121.9,scale); #kg/hr
    CR_Cost=np.multiply(6887460,scale); #USD

    CR_USD_per_kg=CR_Cost*(1/CR_Eff)*(1/CR_LS)*(1/870); #$/kg


    # Hydrogen Redux

    HR_LS=5;  #Years
    HR_Eff=np.multiply(87.64,scale); #kg/hr
    HR_Cost=np.multiply(1141181,scale); #USD

    HR_USD_per_kg=HR_Cost*(1/HR_Eff)*(1/HR_LS)*(1/870); #$/kg

    # Transportation

    #Technology 1

    T1_LS=15;  #Years
    T1_Eff=73.9492; #kg/hr
    T1_Cost=2000000; #USD

    T1_USD_per_kg=T1_Cost*(1/T1_Eff)*(1/T1_LS)*(1/870); #$/kg

    #Technology 2

    T2_LS=10;  #Years
    T2_Eff=24.75; #kg/hr
    T2_Cost=2000000; #USD

    T2_USD_per_kg=T2_Cost*(1/T2_Eff)*(1/T2_LS)*(1/870); #$/kg

    #Technology 3

    T3_LS=5;  #Years
    T3_Eff=31.229; #kg/hr
    T3_Cost=294835.04; #USD

    T3_USD_per_kg=T3_Cost*(1/T3_Eff)*(1/T3_LS)*(1/870); #$/kg

    #Technology 4

    T4_LS=5;  #Years
    T4_Eff=31.229; #kg/hr
    T4_Cost=28520876.19; #USD

    T4_USD_per_kg=T4_Cost*(1/T4_Eff)*(1/T4_LS)*(1/870); #$/kg

    #Technology 5

    T5_LS=5;  #Years
    T5_Eff=31.229; #kg/hr
    T5_Cost=294835040.8; #USD

    T5_USD_per_kg=T5_Cost*(1/T5_Eff)*(1/T5_LS)*(1/870); #$/kg




    ## Products Last edit: Fabrizio

    #Metal: Ti, Si, Mg, Al, Fe, Ca

    #Regolith
    Regolith_Abd=np.multiply(scale,1000);
    Regolith_Cost=25000; #USD/kg

    #MethBucktF,HydrBuckF,MethDrumF,HydrDrumF,MethConvF,HydrConvF
    #Oxygen
    Oxygen_AbdC=[MethBuckF*ProcTime[0][0],MethConvF*ProcTime[1][0],MethDrumF*ProcTime[2][0]];
    Oxygen_AbdH=[HydrBuckF*ProcTime[0][1],HydrConvF*ProcTime[1][1],HydrDrumF*ProcTime[2][1]];
    Oxygen_Cost=11.8; #USD/kg

    #Silicon
    Silicon_Abd= [BuckMetal[1],ConvMetal[1],DrumMetal[1]];
    Silicon_Cost=9.60; #USD/kg

    #Iron
    Iron_Abd= [BuckMetal[4],ConvMetal[4],DrumMetal[4]];
    Iron_Cost=.08; #USD/kg

    #Calcium
    Calcium_Abd= [BuckMetal[5],ConvMetal[5],DrumMetal[5]];
    Calcium_Cost=4.9; #USD/kg

    #Aluminum
    Aluminum_Abd= [BuckMetal[3],ConvMetal[3],DrumMetal[3]];
    Aluminum_Cost=2.25; #USD/kg

    #Magnesium
    Magnesium_Abd= [BuckMetal[2],ConvMetal[2],DrumMetal[2]];
    Magnesium_Cost=3.31; #USD/kg


    ## Revenue and COst Calculations

    Ground_Eff=[[MethBuckF,HydrBuckF],[MethDrumF,HydrDrumF],[MethConvF,HydrConvF]] 
    #Ground_Eff=[MethBuckF,MethDrumF,MethConvF,HydrBuckF,HydrDrumF,HydrConvF] 
    Ground_LS=[BL_LS,BC_LS,DL_LS]
    Product_AbdC=[Regolith_Abd,Oxygen_AbdC,Silicon_Abd,Iron_Abd,Calcium_Abd,Aluminum_Abd,Magnesium_Abd]
    Product_AbdH=[Regolith_Abd,Oxygen_AbdH,Silicon_Abd,Iron_Abd,Calcium_Abd,Aluminum_Abd,Magnesium_Abd]
    Product_Abd=[Product_AbdC,Product_AbdH]
    Product_Cost=[Regolith_Cost,Oxygen_Cost,Silicon_Cost,Iron_Cost,Calcium_Cost,Aluminum_Cost,Magnesium_Cost]
    
    ## Launch Cost Model
    
    #LaunchCost= 67000000; #Falcon 9 Cost
    #GTOMass=8300; #kg mass to GTO
    #MarsMAss= 4020; #kg MAss to Mars
    #LLOMass= (GTOMass + MarsMass)/2 #kg Guess Mass to Moon
    
    #Product_Cost=LaunchCost+LLOMass*Product_Cost*Product_Abd
    
    """
    print("Ground_Eff")
    print(Ground_Eff)
    print("Ground_LS")
    print(Ground_LS)
    print("Product_Abd")
    print(Product_AbdC,Product_AbdH)
    print("Product_Cost")
    print(Product_Cost)
    """

    #Revenue = np.zeros((3,7,2))
    
    Revenue = np.zeros((3,2))
    Value = np.zeros((6,1))


    #print(float(Ground_Eff[0])*float(Ground_LS[0])*(870)*float(Product_Abd[3])*float(Product_Cost[3]))


    for i in range(0,3): #Ground
        for k in range(0,2): #Processing
            TotalVal = 0;
            for j in range(0,7): #Products
                TotalVal += Product_Abd[k][j][i]*Product_Cost[j]
            Revenue[i][k] = Ground_Eff[i][k]*Ground_LS[i]*(870)*TotalVal

    #print("Revenue")
    #print(Revenue)

    #Power Efficiency and Lifespan
    Power_Eff=[GS_Eff,OG30_Eff,OG70_Eff,OS_Eff];
    Power_LS=[GS_LS,OG30_LS,OG70_LS,OS_LS];

    #Cost of infrastructure
    Power_Cost=[GS_USD_per_kJ,OG30_USD_per_kJ,OG70_USD_per_kJ,OS_USD_per_kJ,];
    Ground_Cost=[BL_USD_per_kg,BC_USD_per_kg,DL_USD_per_kg];
    Processing_Cost=[CR_USD_per_kg,HR_USD_per_kg];
    Transportation_Cost=[T1_USD_per_kg,T2_USD_per_kg,T3_USD_per_kg,T4_USD_per_kg,T5_USD_per_kg];
    #Infrastructure_Cost=[Ground][Processing][Power][Transport]
    Infrastructure_Cost= np.zeros((3,2,(4*5)))

    #print(Power_Eff)
    #print(Power_LS)
    #print(Power_Cost)
    #print(Ground_Cost)
    #print(Processing_Cost)
    #print(Transportation_Cost)


    overall_count = 0
    ROICode = []


    for i in range(0,3): #Ground
        for k in range(0,2): #Processing
            counter=0
            for j in range(0,4): #Power
                for h in range(0,5): #Transportation
                    Infrastructure_Cost[i,k,counter] = float(Ground_LS[i])*float(Ground_Eff[i][k])*870*(float(Ground_Cost[i])+float(Processing_Cost[k][i])+float(Transportation_Cost[h]))+float(Power_LS[j])*float(Power_Eff[j])*31536000*float(Power_Cost[j])
                    
                    ROICode.append(str((i,j,h,k)))
                    counter += 1 
                    overall_count += 1


    #print("Infrastructure_Cost")
    #print(Infrastructure_Cost)

    #Only 3 Ground Options

    #Each Individual Ground Option has 7 Revenue Options, and 12 infastructure
    #cost options. Need to Compare each Revenue Option to 12 Infastrcture Cost
    #Options 

    Profit = np.zeros((6*20)) #[[0 for col in range(3*7*12)] for row in range(2)]
    counter = 0
    for i in range(0,3): #Mining Option
        for j in range(0,2): #Processing
            for k in range(0,20): #Cost

                #Revenue_vs_Cost[counter,:]=[Revenue[i,j],Infrastructure_Cost[i,j,k]]; 
                Profit[counter] = Revenue[i,j]-Infrastructure_Cost[i,j,k];
                counter=counter+1;
                #figure(1)
                #plot(Revenue(i,j),Infastructure_Cost(i,k),'o');
                #hold on
                #xlabel('LifeTime Revenue USD')
                #ylabel('Infastructure Cost USD')
                #title('LifeTime Revenue vs Infastrcture Cost')
                
    """            
    print("Revenue v Cost")       
    print(Revenue_vs_Cost)

    Profit= Revenue_vs_Cost[:,0] - Revenue_vs_Cost[:,1]
    print("Profit")
    print(Profit)
    """

    return ROICode, Profit
    
    #return Revenue, Infrastructure_Cost, ROICode, Revenue_vs_Cost, Profit