## Constants  

class TransportAnalysis:
    def  __init__(self):

        self.g = 9.81                   # Earth Gravitational constant [m/s^2] 

        self.gm = 1.625                 # Lunar Gravitational constant [m/s^2] 

        self.lb2kg = 2.20462262         # Divide lb by this to get [kg] 

        self.mph2ms = 2.23693629        # Divide mph by this to get [m/s] 


#_________________________________________________________________________________________________________________ 

#_________________________________________________________________________________________________________________ 

## Saturn V takes Apollo 11/Blue Origin type lander to Moon (Assumes only ascent vehicle is reusable) 
    def ParetoVector1(self,Tt1000,Fuel):

        #print("\n-----------------\n","Tt10000000000000",Tt1000,"\n------------------\n")

        p = 20000                                  # Price to launch one kg to the moon [$] (BASED on Apollo 11) 


        # Apollo 11/Blue Origin type lander 

        m11w = 32500 / self.lb2kg             # Wet mass of apollo 11 lunar lander (w/ propellant and crew) [kg] 

        m11d = 9000 / self.lb2kg              # Dry mass of apollo 11 lunar lander (w/o propellant and crew) [kg] 
         

        mA11w = 4700                     # Wet mass of apollo 11 lunar ascent vechicle (w/ propellant and crew) [kg] 

        mA11d = 2150                     # Dry mass of apollo 11 lunar ascent vechicle (w/o propellant and crew) [kg] 
         

        mD11w = m11w - mA11w             # Wet mass of apollo 11 lunar descent vechicle (w/ propellant and crew) [kg] 

        mD11d = m11d - mA11d             # Dry mass of apollo 11 lunar descent vechicle (w/o propellant and crew) [kg] 
         

        wA11 = self.gm * mA11w                # Gross weight of apollo 11 unar ascent vechicle [N] 

        T11 = 16000                      # Thrust of apollo 11 lunar ascent vechicle [N] 
         

        TtoW11 = T11 / wA11              # Thrust to weight ratio of apollo 11 lunar ascent vechicle 

        minTtoW = 1.5                    # Assumed minimum Thrust to weight ratio of hypothetical lunar ascent vechicle w/ same thrust as apollo 11 one 
         

        MassDiff = (T11 - wA11) / self.gm     # Mass needed to make Thrust to weight ratio of apollo 11 lunar ascent vechicle 0 [kg] 

        # ^ mass of hypothetical lunar ascent vechicle cannot surpass this amount 

        # since we are assume it will have the same thrust as the apollo 11 one 
         

        massAvailable = ((T11/minTtoW)-wA11) / self.gm      # Available mass of lunar ascent vechicle to bring to lunar orbit [kg] 

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

        MassPerHour = SaleOFProductMass2 / (((SaleOFProductMass2+ReusableRequirements2)/Fuel)*1) 

        PowerPerSale = 0 

        
        ParetoVector1 = [MassPerHour, PowerPerSale, ROI] 



        print(ParetoVector1)

        return ParetoVector1
 


#_________________________________________________________________________________________________________________ 

#_________________________________________________________________________________________________________________ 

## Falcon Heavy takes Apollo 11/Blue Origin type lander to Moon (Assumes only ascent vehicle is reusable) 
 
    def ParetoVector2(self,Tt1000,Fuel):
        massAvailableF = 20000                     # FalconH mass to lunar surface [kg] 

        priceFalcon = 100000000                    # Price to launch falcon heavy [$] 

        p = priceFalcon / massAvailableF           # Price to launch one kg to the moon [$] (BASED on Falcon Heavy) 
         

        # Apollo 11/Blue Origin type lander 

        m11w = 32500 / self.lb2kg             # Wet mass of apollo 11 lunar lander (w/ propellant and crew) [kg] 

        m11d = 9000 / self.lb2kg              # Dry mass of apollo 11 lunar lander (w/o propellant and crew) [kg] 
         

        mA11w = 4700                     # Wet mass of apollo 11 lunar ascent vechicle (w/ propellant and crew) [kg] 

        mA11d = 2150                     # Dry mass of apollo 11 lunar ascent vechicle (w/o propellant and crew) [kg] 
         

        mD11w = m11w - mA11w             # Wet mass of apollo 11 lunar descent vechicle (w/ propellant and crew) [kg] 

        mD11d = m11d - mA11d             # Dry mass of apollo 11 lunar descent vechicle (w/o propellant and crew) [kg] 
         

        wA11 = self.gm * mA11w                # Gross weight of apollo 11 unar ascent vechicle [N] 

        T11 = 16000                      # Thrust of apollo 11 lunar ascent vechicle [N] 
         

        TtoW11 = T11 / wA11              # Thrust to weight ratio of apollo 11 lunar ascent vechicle 

        minTtoW = 1.5                    # Assumed minimum Thrust to weight ratio of hypothetical lunar ascent vechicle w/ same thrust as apollo 11 one 
         

        MassDiff = (T11 - wA11) / self.gm     # Mass needed to make Thrust to weight ratio of apollo 11 lunar ascent vechicle 0 [kg] 

        # ^ mass of hypothetical lunar ascent vechicle cannot surpass this amount 

        # since we are assume it will have the same thrust as the apollo 11 one  
         

        massAvailable = ((T11/minTtoW)-wA11) / self.gm      # Available mass of lunar ascent vechicle to bring to lunar orbit [kg] 

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

        MassPerHour = SaleOFProductMass2 / (((SaleOFProductMass2+ReusableRequirements2)/Fuel)*1) 

        PowerPerSale = 0 

        ParetoVector2 = [MassPerHour, PowerPerSale, ROI] 

        print(ParetoVector2) 

        return ParetoVector2

#_________________________________________________________________________________________________________________ 

#_________________________________________________________________________________________________________________ 

## Starship takes Apollo 11/Blue Origin type lander to Moon (Assumes only ascent vehicle is reusable) 

    def ParetoVector3(self,Tt1000,Fuel):

        massAvailableS = 100000                    # Starship mass to lunar surface [kg] 

        priceStarship = 2000000                    # Price to launch falcon heavy [$] 

        p = priceStarship / massAvailableS         # Price to launch one kg to the moon [$] (BASED on Starship)  

        # Apollo 11/Blue Origin type lander 

        m11w = 32500 / self.lb2kg             # Wet mass of apollo 11 lunar lander (w/ propellant and crew) [kg] 

        m11d = 9000 / self.lb2kg              # Dry mass of apollo 11 lunar lander (w/o propellant and crew) [kg] 
         

        mA11w = 4700                     # Wet mass of apollo 11 lunar ascent vechicle (w/ propellant and crew) [kg] 

        mA11d = 2150                     # Dry mass of apollo 11 lunar ascent vechicle (w/o propellant and crew) [kg] 
         

        mD11w = m11w - mA11w             # Wet mass of apollo 11 lunar descent vechicle (w/ propellant and crew) [kg] 

        mD11d = m11d - mA11d             # Dry mass of apollo 11 lunar descent vechicle (w/o propellant and crew) [kg] 
         

        wA11 = self.gm * mA11w                # Gross weight of apollo 11 unar ascent vechicle [N] 

        T11 = 16000                      # Thrust of apollo 11 lunar ascent vechicle [N] 


        TtoW11 = T11 / wA11              # Thrust to weight ratio of apollo 11 lunar ascent vechicle 

        minTtoW = 1.5                    # Assumed minimum Thrust to weight ratio of hypothetical lunar ascent vechicle w/ same thrust as apollo 11 one 


        MassDiff = (T11 - wA11) / self.gm     # Mass needed to make Thrust to weight ratio of apollo 11 lunar ascent vechicle 0 [kg] 

        # ^ mass of hypothetical lunar ascent vechicle cannot surpass this amount 

        # since we are assume it will have the same thrust as the apollo 11 one  
         
        massAvailable = ((T11/minTtoW)-wA11) / self.gm      # Available mass of lunar ascent vechicle to bring to lunar orbit [kg] 

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

        MassPerHour = SaleOFProductMass2 / (((SaleOFProductMass2+ReusableRequirements2)/Fuel)*1) 

        PowerPerSale = 0 

        ParetoVector3 = [MassPerHour, PowerPerSale, ROI] 

        print(ParetoVector3)  

        return ParetoVector3

#_________________________________________________________________________________________________________________ 

#_________________________________________________________________________________________________________________ 

## Starship takes Starship type lander to Moon (Assumes vehicle is fully reusable) 

    def ParetoVector4(self,Tt1000,Fuel):

        massAvailableS = 100000                    # Starship mass to lunar surface [kg] 

        priceStarship = 2000000                    # Price to launch falcon heavy [$] 

        p = priceStarship / massAvailableS         # Price to launch one kg to the moon [$] (BASED on Starship) 
         

        # Starship type lander 

        m11w = 2000000 / p                     # The variable name means nothing here this value is just meant to give cost of the  

                                               # ^starship vehicle for the print statement 
                                                 

        massAvailable = 100000                 # Available mass of lunar ascent/descent vechicle to bring to lunar orbit [kg] 

        propWeightA11 = (1200*1000) / (self.g/self.gm)   # Weight of propellant of Starship vechicle needed per lunar launch [kg] 

        propWeightD11 = 0                      # There is no descent vehicle so set propellant weight to 0 [kg] 
         

        IceBoxCost2a = 0                                           # Extra cost for The Ice Box per sale [$] 

        IceBoxCost2b = m11w * p                                    # Extra ONE TIME cost for The Ice Box [$] 

        CustomerCost2a = 0                                         # Extra cost for customer per sale [$] 

        CustomerCost2b = 0                                         # Extra ONE TIME cost for customer [$] 

        SaleOFProductMass2 = massAvailable                         # Amount of product that can be transported per sale [kg] 

        ReusableRequirements2 = propWeightA11 + propWeightD11      # Amount of fuel needed per sale [kg] 
         

        ROI = IceBoxCost2a + IceBoxCost2b 

        MassPerHour = SaleOFProductMass2 / (((SaleOFProductMass2+ReusableRequirements2)/Fuel)*1) 

        PowerPerSale = 0 

        ParetoVector4 = [MassPerHour, PowerPerSale, ROI] 

        print(ParetoVector4) 

        return ParetoVector4

#_________________________________________________________________________________________________________________ 

#_________________________________________________________________________________________________________________ 

## Starship takes Spin Launch to Moon 
    def ParetoVector5(self,Tt1000,Fuel):

        massAvailableS = 100000                    # Starship mass to lunar surface [kg] 

        priceStarship = 2000000                    # Price to launch falcon heavy [$] 

        p = priceStarship / massAvailableS         # Price to launch one kg to the moon [$] (BASED on Starship) 


        spinw = 100000                                         # Weight of spin launch without vacuum infrastructure [kg] 

        EarthMassLaunch = 200                                  # Mass of object launched on Earth [kg] 

        MoonMassLaunch = (self.g/self.gm) * EarthMassLaunch              # Mass of object launched on Moon [kg] 

        WeightOFequipment = 200                                # Weight of equipment needed to control product after launch [kg] 

        V = 5000 / self.mph2ms                                      # Velocity of projectile at launch [m/s] 

        Energy = (1/2) * MoonMassLaunch * V**2                  # Energy needed per launch [J]  

        # Using Spin Launch 

        IceBoxCost2a = 0                                           # Extra cost for The Ice Box per sale [$] 

        IceBoxCost2b = spinw * p                                   # Extra ONE TIME cost for The Ice Box [$] 

        CustomerCost2a = 0                                         # Extra cost for customer per sale [$] 

        CustomerCost2b = 0                                         # Extra ONE TIME cost for customer [$] 

        SaleOFProductMass2 = MoonMassLaunch - WeightOFequipment    # Amount of product that can be transported per sale [kg] 

        ReusableRequirements2 = Energy                              # Amount of energy needed per sale [J] 


        ROI = IceBoxCost2a + IceBoxCost2b 

        MassPerHour =   Fuel #SaleOFProductMass2 / (((SaleOFProductMass2)/Tt1000)*24) 

        PowerPerSale = ReusableRequirements2 * (3/4)

        ParetoVector5 = [MassPerHour, PowerPerSale, ROI] 

        print(ParetoVector5)

        return ParetoVector5