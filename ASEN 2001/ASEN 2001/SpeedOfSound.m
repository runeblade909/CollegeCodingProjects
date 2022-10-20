%%Speed of Sound (Mach) at given temperature calculator
%Equations from https://www.grc.nasa.gov/www/k-12/airplane/sound.html

function mach = SpeedOfSound(temp)
y = 1.42235; %ratio of specific heats for air

R = 286; %Gas constant for air m^2/s^2/K
temp2 = temp;
w = input('Is this in fahrenheit(0) or celcius(1)?\n');
    if w == 0
    temp = (temp-32) * (5/9);
    temp = temp + 273.15;
    elseif w == 1
    temp = temp + 273.15;
    else
        error('Please type "0" for fahrenheit or "1" for celcius');
    end
    
        
mach = round(sqrt(y*R*temp))/343;
fprintf('The speed of sound at %.2f of given degrees is mach %.2f',temp2,mach); 
end