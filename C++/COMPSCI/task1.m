% rowvec1 = [2, 4, 6, 8, 10]
% rowvec2 = [-10, -20, -30]
% colvec1 = [15; 30; 45; 60]
% colvec2 = [5; 3; 1; -1; -3; -5]

% By Zak Reichenbach
% 9/17/2019
% Make these with colon operator and linspace.
clc
clear

rowvec1col=2:2:10
rowvec1lin= linspace(2,10,5)

rowvec2col=-10:-10:-30
rowvec2lin= linspace(-10,-30,3)

colvec1col=[15:15:60]'
colvec1lin=linspace(15,60,4)'

colvec2col=[5:-2:-5]'
colvec2lin=linspace(5,-5,6)'
