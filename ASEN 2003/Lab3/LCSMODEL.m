function [Beta,theta,velocity] = LCSMODEL(r , d, l, theta, w)
 w = w*((pi)/180);

    Beta = asind((d-r.*sind(theta))./l);
%RADIANS
    velocity = -w.*r.*(sind(theta)+cosd(theta).*tand(Beta));
    



end