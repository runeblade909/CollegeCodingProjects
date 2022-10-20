%largestLast.m
clear
clc

%we want to find the largest value in the code, and move it to the end of
%the vector

v= [9,5,7,1,3,0,0];
N=length(v);

disp(v);
%The outter loop makes the process repeat over and over so that it sorts
%all of the numbers, not just the single largest
for j=1:N
    for i=1:N-j
       if v(i)>v(i+1)
           %swap
           temp = v(i);
           v(i) = v(i+1);
           v(i+1) = temp;
       end
       pause(.5)
       disp(v)
    end
end