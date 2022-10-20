%Task1.m
%By: Zak Reichenbach
%9/26/2019

%Create a matrix with the following values without hard coding the numbers
%directly

clc
clear

A1=13:-3:7;
A2=1:4:9;
A3=5:10:25;

%I know it is expected of me to use the disp() function to display my
%results, but I know if I don't suppress the answers they will pop up
%anyways, so imma do that instead.

A=[A1;A2;A3]

m1=A(3,3)

m2=A(1:3,3)

m3=A(2:3,1:3)



