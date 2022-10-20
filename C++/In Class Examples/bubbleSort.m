%Lab4.m
function [sorted]= bubbleSort(x)
disp(x);

N=length(x);
for j=1:N
    for i=1:N-j
        if x(i)<x(i+1)
        %swap
        temp=x(i);
        x(i)=x(i+1);
        x(i+1)=temp;
        end
        pause(.5)
        disp(x)
    end

end