function mountain = gaussianInitialize(center,E)

syms x y;
x = x - center(1);
y = y - center(2);
r = sqrt(x.^2+y.^2);

mountain = exp(-(E*r).^2);

end