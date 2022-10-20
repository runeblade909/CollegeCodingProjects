%unique_elements

function [diff_prod,r_vec,c_vec] = unique_elements(X,Y)
%here we have 2 matricies, we want to find all of the values that are not
%equivalent, and then find the product of them, and their indicies.
clc


[row,column]=size(X);

logical(X==Y);
    for r=1:row        
         for c = 1:column
            
           if X(c)~=Y(c)
               diff_prod=X(c)*Y(c)
           
            end

         end
    end
end