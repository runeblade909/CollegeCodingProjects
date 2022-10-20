%By Zak Reichenbach
%9/10/2019



function [out] = addFunc(x,y)
%This function will add two numbers
%After the function is given its header, you need to make values for the
%two arguments that you put in the function. I did them in the form of
%inputs from the user so that the function would work properly and be a
%useful function rather than a function stub with dumby data that does'nt
%work.
x=input('Number 1:');
y=input('Plus Number 2:');

out=x+y;
%For now we just assign some dummy data to the output variable, is what I
%would've done if I wanted to follow the exact assignment, but I wanted to
%challenge myself to make code that worked. So I made the out do the
%process that the name of the function promised.
end
