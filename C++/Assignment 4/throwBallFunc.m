function [x,y,hitDistance,hitTime,xt,yt] = throwBallFunc(v,theta)


h = 1.5;

gravity = 9.8;

t = linspace(0,20,10000);

xt = v*cos(theta * (pi/180))*t;

yt = h + v*sin( theta * (pi/180))*t - (.5*gravity*t.^2);

%find the index in the time vector right before the value goes negative
t_hit = find(yt>0& yt<0.005);

x = xt(t_hit);

y = yt(t_hit);

hitDistance = x;

hitTime = t(t_hit);

[x,y,hitDistance,hitTime];
end