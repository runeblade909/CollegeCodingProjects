function [theta_exp, w_exp, v_exp, time] = LCSDATA(filename)
%% Column: Time, Wheel Position(deg), Slide Position(mm), Wheel Speed(deg/s), Slide Speed(mm/s), Actual Sample Time (ms)  

%Load File
file = load(filename);
%Split up file

z= 218; %This is how long we want all of our files to be.
time = file(:,1);
theta_exp = file(:,2);
w_exp = file(:,4);
v_exp = file(:,5)/10; %mm to cm

%Theta Scaling back to zero
theta_exp = theta_exp - 360 * floor(theta_exp(1)/360);

end