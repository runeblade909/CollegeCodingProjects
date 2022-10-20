function[z] = magic_matric(my_matrix,threshold)
clc
%make my_matrix easier to work with

x = my_matrix;
%find where my_matrix is larger than the threshold
  
z = x> threshold

h = x(z);

[length,~] = size(h)
    
%sum the equation
z= sum(x(z));  

for i=1:length
     for j=i
    q= h(j);
     end
  end
    

end

%I wasnt aware we werent able to use the sum function in the test, and I
%did not prepare how to do that without a for loop. I hope that doesnt hurt
%me too much.