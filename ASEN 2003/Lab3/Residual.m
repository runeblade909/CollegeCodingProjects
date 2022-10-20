function [Residuals, Res_mean, Res_std, Resdiuals_2, time_2] = Residual(V_exp,V_theo,time)
% This function calculates the absolute value of the difference between
% the experimental and theoretical velocity values


Residuals = abs(V_exp-V_theo);

Res_mean = mean(Residuals);

Res_std = std(Residuals);

%% Sort out outliers

j = 1;

for i = 1:length(Residuals)

    if Residuals(i) <= Res_mean + Res_std && Residuals(i) >= Res_mean - Res_std 
        
        Resdiuals_2(j) = Residuals(i);
        
        time_2(j) = time(i);
        
        j = j+1;  
        
    end
end


