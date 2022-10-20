function mountain = inversemultiquadraticInitialize(center,E)

syms x y;
x = x - center(1);
y = y - center(2);
r = sqrt(x.^2+y.^2);

mountain = 1./(sqrt(1+(E*r).^2));

end