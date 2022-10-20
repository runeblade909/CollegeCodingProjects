function mountain = inversequadraticInitialize(center,E)

syms x y;
x = x - center(1);
y = y - center(2);
r = sqrt(x.^2+y.^2);

mountain = 1./(1+(E*r).^2);

end