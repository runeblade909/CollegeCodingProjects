%% Properties at Altitude Calculator
%Equations from https://www.grc.nasa.gov/www/k-12/airplane/atmos.html
function [temp,pressure,density] = PropertiesAtAlt(h)
    if h > 82345
        temp = -205.05 + .00164*h;
        pressure = 51.97 * ((temp + 459.7)/389.98)^(-11.388);
    elseif h<82345 && h> 36152
        temp = -70;
        pressure = 473.1 * exp(1.73-.000048*h);
    elseif h<36152
        temp = 59 - .00356*h;
        pressure = 2116 * ((temp + 459.7)/518.6)^(5.256);
    end
     density = pressure /(1718*(temp+459.7));
    fprintf('You now have the temperature in FAHRENHEIT after inserting altitude in FEET :%.2f\n',temp);
    fprintf('You also now have pressure in POUNDS per SQUARE FEET: %.2f\n',pressure);
    fprintf('Pressure in POUNDS per SQUARE INCH: %.2f\n',pressure/144);
    fprintf('Finally, you have density in SLUGS per CUBIC FEET: %f\n',density);
    

end
