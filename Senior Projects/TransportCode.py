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