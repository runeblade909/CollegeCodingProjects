function [Residuals, Res_mean, Res_std, Residuals_2, time_2, Res2_mean, Res2_std] = Residual(V_exp,V_theo,time)
% Purpose: calculates the absolute value of the difference between
% the experimental and theoretical velocity values
% Author: Jacob Starkel
% Date Completed: 3/12/2021

Residuals = V_exp-V_theo;

Res_mean = mean(Residuals);

Res_std = std(Residuals);

%% Sort out outliers

j = 1;

for i = 1:length(Residuals)

    if Residuals(i) <= Res_mean + 1.5*Res_std && Residuals(i) >= Res_mean - 1.5*Res_std 
        
        Residuals_2(j) = Residuals(i);
        
        time_2(j) = time(i);
        
        j = j+1;  
        
    end
end

Res2_mean = mean(Residuals_2);

Res2_std = std(Residuals_2);
