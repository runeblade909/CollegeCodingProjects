function [e_time] = E_t(BounceTime)
%This calculates the coefficient of resititution with time values

    for i = 1:length(BounceTime)-1
    e_time(i) = BounceTime(i+1)/BounceTime(i);
    end
    
end