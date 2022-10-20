function [E_Th] = E_th(BounceTime,BounceHeight)
g = 9.81; %[m/s^2]
Ts = sum(BounceTime);
h0 = BounceHeight(1);
E_Th = (Ts - sqrt((2*h0)/g))/(Ts + sqrt((2*h0)/g));

end