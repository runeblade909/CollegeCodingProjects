function mountain = gaussian(centerX,centerY,E)

[x,y] = meshgrid(0:0.5:5,0:5);

mountain = exp(E*(sqrt((x-centerX)^2+(y-centerY)^2)))^2;


end