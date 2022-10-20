function [e_height] = E_h(BounceHeight)
%This calculates the coefficient of resititution with height values

    for i = 1:length(BounceHeight)-1
    e_height(i) = sqrt(BounceHeight(i+1)/BounceHeight(i));
    end
    
end