% In a script task2.m, generate the following matrix M using two vectors vec1 = [1, 2, 4] and vec2 = [5,
% 6, 8] for the first and second rows and the colon operator for the third row.
% M =
% 1 2 4 5 6 8
% 5 6 8 1 2 4
% 1 3 5 7 9 11

% By Zak Reichenbach
% 9/17/2019
% Do as it says

clc
clear


vec1=[1,2,4];
vec2=[5,6,8];

x=[vec1,vec2;vec2,vec1;linspace(1,11,6)]