function mountain = multiquadraticInitialize(center,E)

syms x y;
x = x - center(1);
y = y - center(2);
r = sqrt(x.^2+y.^2);

mountain = sqrt(1+(E*r).^2);

end