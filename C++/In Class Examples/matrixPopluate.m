%matrixPopluate.m
clear
clc

N=10;

m=zeros(N);

for j = 1:N
    for i = 1:N
        m(j,i) = i*j;
        image(m)
        pause(.1)
    end
end
%image(m)