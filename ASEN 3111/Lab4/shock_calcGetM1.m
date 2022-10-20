function [ M1 ] = shock_calcGetM1( p2op1 )

g=1.4;
M1 = sqrt((p2op1 -1)^((g+1)/(2*g))+1);

end