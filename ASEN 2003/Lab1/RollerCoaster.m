%% House keeping

%unit circule starts on left, goes negative in counterclockwise direction
%for positive angles

   clc
   clear
   close all
   
%% Constants
h0 = 125;%m
g = 9.8; %m/s^2
m = 700; %kg
   
%% Calculate velocity of particle for any height
v = zeros(125,1);
h = 1:125;
 for i = 1:125
     v(i) = sqrt(2*g*(h0-h(i)));
 end
 v=v';
%% Intial Drop (For the first quarter circle)
h1 = h0-20; %Drop 20 meters for speed
v1 = v(h1);
drop = h0:-1:h1;

%% Quarter Circle
r1 = 15; % [m]
C1 = [r1 h1] ; %Center of circle 
th1 = 0:1:90 ; %Change in theta

%Circle Equations
xC1 = C1(1)- r1*cosd(th1) ; 
zC1 = C1(2)- r1*sind(th1) ;
%Ending Coordinates
h2 = zC1(end);
x2 = xC1(end);

%pi/2 = 90 degrees  circumference is 2*pi*r
length1 = (pi/2)*r1;

%G-Forces
g1_max = g_Circle(v(h2), th1, g, r1);
g1_min = g_Circle(v(h1),th1,g,r1);


%% Linear Region 1
xline = x2:1:x2+10; %10 m long
yline = zeros(1,length(xline)); %constant y plane
zline = h2*ones(1,length(xline)); % constant z plane



%% Loop 1
r2 = 15; %[m]
C2 = [xline(end) yline(end) zline(end)+r1]; %center
th2 = 90:4:450; %special angles for whole circle

xC2 = C2(1)-r2*cosd(th2); 
zC2 = C2(3)-r2*sind(th2); 
yC2 = zeros(1,length(xC2));

%End of second circle
h4 = zC2(end);
x4 = xC2(end);
y4 = yC2(end);

%length full circle
length2 = 2*pi*r2;

%G-Forces
g2_max = g_Circle(v(zC2(1)), th2, g, r2);
g2_min = g_Circle(v(max(zC2)),th2,g,r2);

%% Linear Region 1
line1=(x4:1:x4+50); % 50m long

x4 = line1(end);

g3 = 1*ones(1,length(line1));

%% Banked turn (Loop2)
bank_angle = 45;
r3 = v(x4(end))^2/g*cot(bank_angle);  % Radius needed to have 0 Lateral force on the banked turn
C3 =[x4 y4-r3 h4]; %on horizontal plane, center of banked turn circle
th3 = 180:2:360;
%Line Equations
xC3 = C3(1) - r3*sind(th3); 
yC3 = C3(2) - r3*cosd(th3); 
zC3 = h4*ones(1,length(yC3));
%End Coordinates
x5 = xC3(end);
y5 = yC3(end);
h5 = zC3(end);


%length 
length3 = pi*r3;
%G_Forces
g4 = g_Bank(bank_angle)*ones(round(length3),1);

%% Next Linear Region

line2 = (x5:-1:x5-15); %15m long

g5 = 1*ones(length(line2),1);

x6 = line2(end);

%% Zero G Drop 

   %Setting Inital Conditions for projectile motion after bank
   theta = 0;  %Flat launch angle from banked turn
   th_Parabola = -90:-(91/51):-180;
   th_Parabola = flip(th_Parabola);
   v0 = v(h5); %Starting velocity of banked turn
   g = 9.81;   
   h_change = 50; %Free fall trajectory that falls 50m
   h_vec = 0:h_change;
   
   [x0g, t_change] = zerog(theta,v0,g,h_change);
   
   t_change = t_change';
   
   h = x0g(end);
   k = h_change;
   
   a = k/(h^2);
   
   %Function for projectile motion
   y = @(x) -a*((x-h)^2) + k;
   %Plugging in for projectile motion values
   y0g = -a*((x0g-h).^2)+ k;
   x7 = x0g+(x6-x0g(end));
   y7 = y5*ones(1,length(y0g));
   h7 = (y0g+(40))';
   
   syms x
   y = -a*((x-h).^2)+k;
%derivative of zero g function
   d_y = -2*a*(x-h);
%2nd derivative of zero g function
   d_y2 = diff(d_y);
%derivative at the end of the drop
   dy_ds = sym2poly(subs(d_y,x,x0g(1)));
   
   dz_ds = subs(d_y,x,x0g);
  
   for i = 1:length(dz_ds)
   dz1_ds(i,1) = sym2poly(dz_ds(i,1));
   end
   dy2_ds = sym2poly(d_y2);
   
   dz_ds = dz_ds';
   
   dz1_ds=dz1_ds';
   
%  r_Parabola = (1+(dz1_ds).^2).^(3/2)/((-(dy2_ds)));

   r_Parabola = sqrt((x0g-x0g(end)).^2+(y0g).^2);
   r_Parabola = r_Parabola';
   
 v_Parabola = v(round(h7));
 g6 = g_Parabola(v_Parabola,th_Parabola,g,r_Parabola); 
 
 
 w = sqrt(1+(-2*a*(x-h))^2);
 
 %Length of parabola
 length4 = sym2poly(int(w,0,x0g(end)));
 
 %G-Force
 
 
 
 
%    figure(2)
%    fplot(y);
%    ylim([0 60]);
%    xlim([-10,30]);
%    xlabel('Horizontal Distance');
%    ylabel('Vertical Distance');
%    title('Zero G Parabola');

%% Ramp
% NUMBER FOR HEIGHT OF LAST SEMI CIRCLE IS 11.6978, x distance needed to be traveled
% to drop to zero
t =28.2862/dy_ds;

x_ramp = 0:.05:t;

%deltas * deltay/deltas = deltay
y_ramp = x_ramp*dy_ds;

x_ramp_adjusted = x_ramp +(x7(1)-x_ramp(end));

y_ramp_adjusted = y_ramp +(h7(1)-y_ramp(end));

x8 = x_ramp_adjusted(1);
h8 = y_ramp_adjusted(1);

%Length
length5 = sqrt(x_ramp(end)^2+y_ramp(end)^2);

%G-Force
g7 = cos(atan(dy_ds/1))*ones(round(length5),1);


%% Last exit circle

r4 = 50;
C4 = [x8-r4 y7(end) h8]; %Center
th4 = 90:1:130; %Angle needed to be covered, found experimentally
%Equations
xC4 = C4(1) - r4*cosd(th4);
yC4 = y5*ones(1,length(xC4));
zC4 = C4(3) - r4*sind(th4);
%Analysis
height = zC4(end) - zC4(1); %The height of the quarter circle
height2 = h8 - y_ramp_adjusted(end); %The height needed to be covered by the ramp to reach zero
%End points w/ shift
x9 = xC4(1);
y9 = yC4(1);
z9 = zC4(1);
z9_adjusted = zC4(1)+(h8(1)-zC4(end));

length6 = (length(th4)*(pi/180))*r4;


g8_max = g_Circle(v(round(zC4(1)+(h8(1)-zC4(end)))+1), th4, g, r4);
g8_min = g_Circle(v(round(zC4(end)+(h8(1)-zC4(end)))+1), th4, g, r4);

%% Breaking region

xline2 = x9+(x8(1)-xC4(end)):-1:x9-200+(x8(1)-xC4(end));
%Shifted line that is 50m long to enable a good braking window

vbreak = v(round(z9_adjusted*ones(1,length(xline2)))+1) - linspace(0,49.2991,length(xline2));


%% Total Length

total_length = length(drop)+length1+length(xline)+length2+length(line1)+length3+length(line2)+length4+length5+length(xline2);

%% Gs plot

figure(2)
hold on
plot(th1,g1_max)
plot(th1,g1_min)
title('Gs of quarter Circle')
xlabel('Angle [Degrees]')
hold off
figure(3)
hold on
plot(th2,g2_max)
plot(th2,g2_min)
title('Gs of Loop')
xlabel('Angle [Degrees]')
hold off
figure(4)
plot(1:length(g4),g4)
title('Gs of Banked Turn')
xlabel('Length of bank [m]')
figure(5)
plot(1:round(length5),g7*ones(1,round(length5)))
title('Gs of Ramp')
xlabel('Length of ramp [m]')
figure(6)
hold on
plot(th4,g8_max)
plot(th4,g8_min)
xlabel('Angle [Degrees]')
title('Gs of Exit Circle')
hold off

%% Plot
figure(1)

scatter3(zeros(1,length(drop)),zeros(1,length(drop)),drop,1,v(drop))
hold on
scatter3(xC1,zeros(1,length(zC1)),zC1,1,v(round(zC1)));
scatter3(xline,yline,zline,1,v(zline));
scatter3(xC2,yC2,zC2,1,v(round(zC2)));
scatter3(line1,(y4)*ones(1,length(line1)),h4*ones(1,length(line1)),1,v(h4*ones(1,length(line1))))
scatter3(xC3,yC3,zC3,1,v(zC3));
scatter3(line2,y5*ones(1,length(line2)),h5*ones(1,length(line2)),1,v(h5*ones(1,length(line2))))
scatter3(x7 ,y7, h7,1,v(round(h7)));
scatter3(x_ramp_adjusted,y7(1)*ones(1,length(x_ramp)),y_ramp_adjusted,1,v(round(y_ramp_adjusted)))
scatter3(xC4+(x8(1)-xC4(end)),yC4,zC4+(h8(1)-zC4(end)),1,v(round(zC4+(h8(1)-zC4(end)))+1))
scatter3(xline2,y7(1)*ones(1,length(xline2)),z9_adjusted*ones(1,length(xline2)),1,vbreak)
q = colorbar;
ylabel(q, 'Velocity [m/s]')

title('Our "Safe" Rollercoaster')
xlabel('x[m]')
ylabel('y[m]')
zlabel('h[m]')
hold off



function [x1, t_change] = zerog(theta,v0,g,h_change)
% Preallocate vectors for change in x and change in t
t_change = zeros(h_change+1,1);
x1 = zeros(h_change+1,1);

h_vec = 0:h_change;

% First calculate change in time at each y increment
    for i = 1:h_change+1
        t_change(i) = ((-1)*v0*sind(theta) + sqrt((v0*sind(theta))^2 - 4*(.5*g)*(-1)*(h_vec(i))))/(2*.5*g);
    end

% Calculate change in x at each change in y
    for j = 1:h_change+1
        x1(j) = v0*cosd(theta)*t_change(j); 
    end

end

function [gs] = g_Circle(velocity, theta, gravity, radius)

gs = (velocity^2)/(gravity*radius) - sin(theta);

end

function [gs] = g_Bank(theta)
%we want zero g's in the lateral force, so the term simplifies drastically
%for the bank

gs = 1/cos(theta);

end

function [gs] = g_Parabola(velocity, theta, gravity, radius)

gs = sin(theta) + ((velocity.^2)./(radius.*gravity));

end
