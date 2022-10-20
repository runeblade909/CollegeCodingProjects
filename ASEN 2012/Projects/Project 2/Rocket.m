%% ASEN 2012 Project 2: Bottle Rocket Modeling
%Zak Reichenbach
%11/11/2020

%% House keeping
clc;clear;

%% Governing systems

%Newtons second law
%%Sum(Forces) = M*a = M*[Ax;Az] = M*Vdot = Thrust - Drag + M*g (1)

%%Drag = q*Cd*Ab = (1/2*roe*V^2)*Cd*Ab                         (2)
%Cd is the drag coefficient
%Ab is the cross-sectional area of the front of the bottle

%% Water
%%(Ptot/Pair) = (VolAir/VolTot)^1.4                            (3) 

%%MassFlow = Cd*RoeW*At*Ve                                     (4)
%Cd is discharge Coefficient
%RoeW is density of water
%At is the area of the throat
%Ve is the velocity of exhaust

%%Thrust = MassFlow*Ve +(Pe - Pa)*At                           (5)
%Pa is ambient pressure

%Using Bernoulli's for velocity

%%Ve = sqrt(2(Pair - Pa)/RoeW)                                (6-7)

%Since Pe = Pa, second term of (5) equals 0.

%%Thrust = MassFlow(Ve) = 2*Cd*At*(Pair-Pa)                    (8)

%Air pressure Pair decreases with time as the air volume expands, providing
%less thrust. But the air pressure can be computed from equation 3 if
%volume is known.

%%(dV_dt) = Cd*At*Ve = Cd*At*sqrt(2(Pair - Pa)/RoeW) =Cd*At*Ve =
% Cd*At*sqrt(2/RoeW*(Pair(VolAir/VolTot)^1.4-Pa))              (9)

%We use ODE45 here with the initial conditions v = Vair to Vbottle at time t=0.

%%MassFlowofRocket = -MassFlow = -(Cd*RoeW*At*Ve) = -Cd*At*sqrt(2*RoeW*(Pair-Pa)) 
                                                             % (10)  
%%MassRocket = MassBottle + MassWater + MassAir                (11)

%Now for equation (11) in terms of density 

%%MassRocket = MassBottle +RoeW*(VolTot - VolAir) + ((Pair*VolAir)/(R*Tair))
%% Constants

Values = num2cell([9.81;0.8;0.961;0.002;12.1;1.4;1000;2.1;10.5;287;0.15;0.5;50;0.001;300;0.0;45;0.0;0.25;0.5]);
Names =["g","Cd","RoeAmb","VolBot","PAmb","Gamma","RoeWat",'Dt','Db','R','MBot','CD','Pgage','VolWat0','Tair0','V0','Theta','X0','Z0','L']';
C = cell2struct(Values,Names,1);  %Create Struct of Constants

C.Pgage = C.Pgage*6894.76;      %Convert to pascals
C.PAmb = C.PAmb *6894.76;       %Convert to pascals

P0 = C.Pgage + C.PAmb;          %Pabs = Pgage + Patm

At = pi*((C.Dt*10^-2)/2)^2;             %m^2 Area of throat
Ab = pi*((C.Db*10^-2)/2)^2;             %m^2 Area if bottle

VolAir0 = C.VolBot - C.VolWat0; %m^3 Initial Volume of Air

tspan = [0 5];                  %Span of whole test

%Find Initial Mass Values

M0_Water = C.RoeWat*(C.VolWat0);           %(12)
M0_Air = (P0*VolAir0)/(C.R*C.Tair0);       %(12)
M0 = C.MBot + M0_Water + M0_Air;           %(12)

C.Theta = deg2rad(C.Theta);         %Theta ---> Radians
%Initial State
Properties0 = [C.Z0 C.X0 C.V0 C.Theta M0 M0_Air VolAir0];

%% ODE CALCULATION
terminalCond = odeset('Events', @hitGround);

[t,Properties] = ode45(@(t,s) WaterProp(t,s,C,At,Ab,P0,VolAir0,M0_Air), tspan, Properties0, terminalCond);


%% Further Calculations
n = length(Properties);

for i = 1:n
%      P = roe*R*T
%         V_e = sqrt(C.Gamma*C.R*T_e);
%         V_e = M_e*sqrt(C.Gamma*C.R*T_e);   
%          roe_e = p_crit/(C.R*T_e);
%         
%     thrust = ((-C.Cd)*roe_e*AreaThroat*V_e*V_e)+(C.PAmb-p_e)*AreaThroat

                                 
   %% Use volume to calculate pressure for this 
    
end


%% Plot
figure(1)
plot(Properties(:,2),Properties(:,1));
title('Rocket Propulsion Trajectory')
xlabel('Distance [m]')
ylabel('Height [m]')
figure(2)
plot(t,Properties(:,3));
title('Velocity over Time')
xlabel('Time [S]')
ylabel('Velocity [m/s]')
figure(3)
plot(t,Properties(:,7));
title('Air Volume')
xlabel('Time [S]')
ylabel('Volume [m^3]')
xlim ([0 0.25])
ylim ([1*10^-3 2*10^-3])



%% Stand parameters
StandLength = C.L; %Meters
StandHeight = C.Z0;
StandXPos = C.X0;


%% Functions 
 function state = WaterProp(t,state_i,C,AreaThroat,AreaBottle,P0,Vol0,M0_Air)
 
 
 %State Vector
 state = zeros(7,1);
 Z = state_i(1);                          %Z position [m]
 X = state_i(2);                          %X position [m]
 V = state_i(3);                          %Velocity [m/s]
 Theta = state_i(4);                      %Theta [radians]
 M_Tot = state_i(5);                      %Mass of Rocket [kg]
 M_Air = state_i(6);                      %Mass of Air    [kg]
 VolAir = state_i(7);                     %Volume of Air  [m^3]


 
    %Water Propulsion 
    if(VolAir < C.VolBot)
       %Change in Vol
       dv_dt = C.Cd*AreaThroat*sqrt((2/C.RoeWat)*(P0*((Vol0/VolAir)^C.Gamma)-C.PAmb));
       %Change in pressure
       AirPressure = P0*(Vol0/VolAir)^C.Gamma;
       %Change in thrust
       Thrust = 2*C.Cd*AreaThroat*(AirPressure-C.PAmb);
       %Mass rate of change of total mass
       M_dot = (-C.Cd)*C.RoeWat*AreaThroat*sqrt((2*(AirPressure-C.PAmb))/C.RoeWat);
       %Mass rate of change of air
       M_dot_Air = 0;
    else
        %Pressure of air when all water is gone
        p_end = P0*(Vol0/C.VolBot)^C.Gamma;             %(13)
        %Temperature of air when all water is gone
        t_end = C.Tair0*(Vol0/C.VolBot)^(C.Gamma-1);    %(13)
        %Pressure at end
        AirPressure = p_end*(M_Air/M0_Air)^C.Gamma;     %(14)
        
    end

    
    
    %Pressurized Air
    if(VolAir >= C.VolBot)&&(AirPressure>C.PAmb)
        %Air Pressure for pressurized air
        AirPressure = p_end*(M_Air/M0_Air)^C.Gamma;     %(14)
        %Change in vol
        dv_dt = 0;                                      %No more volume Change
        %Critical Pressure Calculation
        p_crit = AirPressure*(2/(C.Gamma+1))^(C.Gamma/(C.Gamma-1)); %(16)
        %Air Density Calculation
        roe_air = M_Air/C.VolBot;                       %(15)
        %Air Temperature Calculation
        T = AirPressure/(roe_air*C.R);                  %(15)
        
        %Flow
        if(p_crit > C.PAmb) %Choked
            %Exit Temp
            T_e = (2/(C.Gamma+1))*T                    %(18)
            %Exit velocity
            V_e = sqrt(C.Gamma*C.R*T_e);                %(17)
            %Exit Density
            roe_e = p_crit/(C.R*T_e);                   %(18)
            %Mass rate of change air
            M_dot_Air = -C.Cd*roe_e*AreaThroat*V_e;     %(23)
            %Mass rate of change total (The only changes are in air)
            M_dot = M_dot_Air;                          
            %Exit Pressure
            p_e = p_crit;
            
        elseif (p_crit <= C.PAmb) %Non-choked           
            %Exit Mach Number
            M_e = sqrt((((AirPressure/C.PAmb)^((C.Gamma-1)/C.Gamma))-1)/((C.Gamma-1)/2));
            %Air Pressure for non-chocked flow
            AirPressure = C.PAmb*(1+((C.Gamma-1)/2)*M_e^2)^(C.Gamma/(C.Gamma-1)); %(19)
            %Exit Temperature
            T_e = T/(1+(((C.Gamma-1)/2)*(M_e^2)));      %(20)
            %Exit Velocity
            V_e = M_e*sqrt(C.Gamma*C.R*T_e);            %(21)
            %Exit Density
            roe_e = C.PAmb/(C.R*T_e);                   %(20)           
            %Mass rate of change air
            M_dot_Air = (-C.Cd)*roe_e*AreaThroat*V_e;     %(23)
            %Mass rate of change total
            M_dot = M_dot_Air;                          %(24)
            p_e = C.PAmb;                                       %How can this assumption be made? no more propulsion?
        end
        %Thrust Calculation
        Thrust = (M_dot_Air*V_e)+(C.PAmb-p_e)*AreaThroat;       %I should be able to replace the first few terms with M_dot_Air?
    end
    %Ballistic Phase
    if(VolAir>= C.VolBot) && (AirPressure < C.PAmb)
        %Empty rocket
        Thrust = 0;
        M_dot = 0;
        M_dot_Air = 0;
        dv_dt = 0;
         
        
    end
 
 
 %Drag Calculation
 
q = (1/2)*C.RoeAmb*V^2;  
Drag = q*C.CD*AreaBottle;


%Velocity Calculation
V_Z = V*sin(Theta);
V_X = V*cos(Theta);

%Acceleration Calculation
% Thrust = 2*C.Cd*AreaThroat*(AirPressure-C.PAmb);
% M_dot = (-C.Cd)*AreaThroat*sqrt(2*C.RoeWat*(AirPressure-C.PAmb));
%  M_dot_Air = 0;
%    dv_dt = C.Cd*AreaThroat*sqrt((2*(((Vol0/VolAir)^C.Gamma*P0)-C.PAmb))/C.RoeWat);
%  
% FThrust = Thrust;
Acceleration = (Thrust - Drag + (M_Tot*C.g))/M_Tot;     % This is very weird when I change the sign of rocket weight
Accel_Z = Acceleration* cos(Theta);
Accel_X = Acceleration* sin(Theta);


%Calculate Theta





if (Z < (C.L*sin(45)))
    theta_dot = 0;
else
    theta_dot =-.5;
end


%State Output Vector
state(1) = V_Z;
state(2) = V_X;
state(3) = Acceleration;
state(4) = theta_dot;
state(5) = M_dot;
state(6) = M_dot_Air;
state(7) = dv_dt;

 
%  
%  Velocity = sqrt(Vx^2+Vz^2);        %How do I account for this?
%  
%  q = (1/2)*C.RoeAmb*Velocity^2;        %What is V?
%  
%  VolAir = C.VolBot -VolWat;      %Vol of water changes, VolWater = VolWat0 - dvdt?
%  
%  Drag = q*C.CD*Ab;
%  
%  AirPressure = (VolAir0/VolAir)^C.Gamma*C.Pgage;   %Is this Pgage for Pair(i)? or is it Pair(i) = Pgage + PAmb
% 
%  MassAir = (AirPressure*VolAir)/(C.R*C.Tair0); %Constant for first phase.
%  
%  MassWat = C.RoeWat*(C.VolBot-VolAir);
%  
%  Thrust = 2*C.Cd*At*(AirPressure-C.PAmb);
%  
%  MassRocket = C.MBot + MassWat + MassAir;
%  
%  Force = Thrust - Drag +MassRocket*C.g;     %Am I supposed to solve this for the acceleration 
%                                             %and integrate for the velocity which can be broken down
%                                             %into Vz and Vx?
%  Ve = sqrt((2*(AirPressure-C.PAmb))/C.RoeWat);
%  
%  dvdt = C.Cd*At*Ve;
 

 
 

 end
 
 
 
 function [value, isterminal, direction] = hitGround(t,state)
 value      =(state(1) < 0);
 isterminal = 1;
 direction =  0;
 end