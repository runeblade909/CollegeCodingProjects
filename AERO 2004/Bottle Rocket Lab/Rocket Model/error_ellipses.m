clc; clear all; close all; 
 
%% Replace this section of code with your real data
% Simulate and plot 100 data points, (YOU SHOULD USE REAL DATA HERE!)

load('landingPoints.mat')

N = simulations; % Number of points to simulate
x = Xf; % Data From ISP Rocket Model
y = Yf;
%%

figure; plot(x,y,'k.','markersize',6)
axis equal; grid on; xlabel('X Down Range [m]'); ylabel('Y Cross Range [m]'); title('Landing Ellipses'); hold on;
 
% Calculate covariance matrix
P = cov(x,y)
mean_x = mean(x);
mean_y = mean(y);
 
% Calculate the define the error ellipses
n=N; % Number of points around ellipse
p=0:pi/n:2*pi; % angles around a circle
 
[eigvec,eigval] = eig(P); % Compute eigen-stuff
xy_vect = [cos(p'),sin(p')] * sqrt(eigval) * eigvec'; % Transformation
x_vect = xy_vect(:,1);
y_vect = xy_vect(:,2);
 
% Plot the error ellipses overlaid on the same figure
plot(1*x_vect+mean_x, 1*y_vect+mean_y, 'b')
plot(2*x_vect+mean_x, 2*y_vect+mean_y, 'g')
plot(3*x_vect+mean_x, 3*y_vect+mean_y, 'r')
