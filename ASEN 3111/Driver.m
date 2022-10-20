%% Zak Reichenbach
%With help from Bennet Grow and Ketan Kanat

clc
clear
close all

load Cp.mat


%% Problem 1 Spinning Sphere
syms R V theta

tao = 2*pi*R*V;

% Simplifies Greatly
CpTheta = 1 - (4*sin(theta)^2 + (2*tao*sin(theta))/(pi*R*V)+(tao/(2*pi*R*V))^2);

cl = -(1/2)*CpTheta*sin(theta);
cd = -(1/2)*CpTheta*cos(theta);


clActual = -(1/2)*int(CpTheta*sin(theta),theta,0,2*pi);
cdActual = -(1/2)*int(CpTheta*cos(theta),theta,0,2*pi);

fprintf('Sectional Coefficient of Lift: %f (2pi)\n',clActual);
fprintf('Sectional Coefficient of Drag: %f\n',cdActual);

for N = 1:10
    %Trap Lift        
clTrap = CompTrapInt(N,cl);
    %Trap Drag
cdTrap = CompTrapInt(N,cd);
    %Simpson Lift    
clSimp = CompSimpInt(N,cl);
    %Simpson Drag
cdSimp = CompSimpInt(N,cd);

clTrapTot(N) = sum(clTrap);
% if clTrapTot(N) > (clActual-(clActual/500)) && clTrapTot(N) < (clActual+(clActual/500))
%    TrapCountCl = N; 
% end
cdTrapTot(N) = sum(cdTrap);
% if(cdTrapTot > cdActual-(cdActual/500) && cdTrapTot < cdActual+(cdActual/500))
%    TrapCountCd = N; 
% end

clSimpTot(N) = sum(clSimp);
% if(clSimpTot > clActual-(clActual/500) && clSimpTot < clActual+(clActual/500))
%    SimpCountCl = N; 
% end
cdSimpTot(N) = sum(cdSimp);
% if(cdSimpTot > cdActual-(cdActual/500) && cdSimpTot < cdActual+(cdActual/500))
%    SimpCountCd = N; 
% end
   
fprintf('I am on interation: %i \n', N)
fprintf('Sectional Coefficient of Lift Trapezoidal: %f (2pi)\n',clTrapTot(N));
fprintf('Sectional Coefficient of Drag Trapezoidal: %f\n',cdTrapTot(N));

fprintf('Sectional Coefficient of Lift Simpson: %f (2pi)\n',clSimpTot(N));
fprintf('Sectional Coefficient of Drag Simpson: %f\n',cdSimpTot(N));


end

%TO DO 
%Coeff of Lift and drag withing 2/10 percent relative error for trap and
%composite


%% Plots

fprintf('Plotting Graphs.....')

%Problem 2 Plot for Coeff. of Pressure
figure()
hold on
grid on
plot(Cp_lower.knots(7:end),-Cp_lower.coefs(:));
plot(Cp_upper.knots(7:end),-Cp_upper.coefs(:));

legend('Lower Wing', 'Upper Wing')
title('NACA 0012 Coeff. of Pressure Distribution')
xlabel('Chord Length X/C')
ylabel('Coeff. of Pressure')

%Trapezoidal Plots

figure()
hold on
grid on
plot(1:N,cdTrap)
plot(1:N,clTrap)
title('Trapezoidal Integration\n')
legend('Coeff. of Drag','Coeff. of Lift')

%Error Trap

figure()
hold on
grid on
plot(1:N,cdTrapTot)
plot(1:N,clTrapTot)
title('Trapezoidal Integration Values over Iterations')

%Simpson Plots

figure()
hold on
grid on
plot(1:N,cdSimp)
plot(1:N,clSimp)
title('Simpson Integration')
%legend('Coeff. of Drag','Coeff. of Lift')

%Error Simp

figure()
hold on
grid on
plot(1:N,cdSimpTot)
plot(1:N,clSimpTot)
title('Simpson Integration Values over Iterations')
%legend('Coeff. of Drag','Coeff. of Lift')

%% Please Remove Me
%close all;

%% Airfoil 

 C = 1.5; %[m]
 a = 9; %[deg]
 V = 40; %[m/s]
 rho = 1.225; %[kg/m^3]
 Pinf = 101.3*10^3; %[pa]
 
 Dynamic =(1/2)*rho*V^2;

 
 T = .01:(1.5/N):1.51;
 
  %Eq for airfoil curvature
 T = .01:(1.5/N):1.51;
 y =(12/100)/.2*C*(0.2969*(T/C).^(1/2)-0.1260*(T/C)-0.3516*(T/C).^2+0.2843*(T/C).^3-0.1036*(T/C).^4);
 syms X;
 ysyms = (12/100)/.2*C*(0.2969*(X/C).^(1/2)-0.1260*(X/C)-0.3516*(X/C).^2+0.2843*(X/C).^3-0.1036*(X/C).^4);
 v = diff(ysyms);
 
 

c = 1.5; 
alpha = 9; 
V_inf = 40;
rho_inf = 1.225; 
p_inf = 101.3*10^3;

%% Lift and Drag per unit span
% determine values using trapezoidal rule

% Find thickness given NACA airfoil name
t = 12/100;

% Input values into trapezoidal rule calculator to find sectional
% coefficients
[Cl_NACA,Cd_NACA,N_NACA] = trap_Problem2(c,alpha,t,Cp_upper,Cp_lower);

% Use sectional coefficients to find forces per unit span
L = .5*Cl_NACA*rho_inf*V_inf^2*c; 
D = .5*Cd_NACA*rho_inf*V_inf^2*c; 

% Obtain accurate value for both by using the value at maximum panel size
L_actual = L(end);
D_actual = D(end);

% Plot results to confirm results are as expected
figure()
hold on
grid on
plot(N_NACA,L)
title('Lift per Unit Span for Each Panel Size Integration')
xlabel('Panel Size (N)')
ylabel('Lift per Unit span (N/m)')
yline(L(end),'--k')
legend('Lift per Unit Span vs. Panel Size','Assumed Solution (at Max Panel Size)','Location','best')
hold off

figure()
hold on
grid on
plot(N_NACA,D)
title('Drag per Unit Span for Each Panel Size Integration')
xlabel('Panel Size (N)')
ylabel('Drag per Unit span (N/m)')
yline(D(end),'--k')
legend('Drag per Unit Span vs. Panel Size','Assumed Solution (at Max Panel Size)','Location','best')
hold off

% Print results
fprintf('The Lift and Drag per unit span were found to be %.4f and %.4f N/m, respectively.\n',L_actual,D_actual)

%% Relative Error
% Find integration points required for relative error to be 5, 1, and 1/10%

% Calculate relative error for every panel size
L_error = abs((L - L_actual)./L_actual);
D_error = abs((D - D_actual)./D_actual);

% Relative error limit
error_5 = .05;
error_1 = .01;
error_10th = .002;

% Find Panel limit for each error
ind_5 = find(L_error <= error_5,1);
N_lim_5 = N_NACA(ind_5);

ind_1 = find(L_error <= error_1,1);
N_lim_1 = N_NACA(ind_1);

ind_tenth = find(L_error <= error_10th,1);
N_lim_tenth = N_NACA(ind_tenth);

fprintf('The number of integration points for 5%%, 1%% and .1%% relative error are %d, %d, and %d.\n',N_lim_5,N_lim_1,N_lim_tenth)


 %Plot airfoil
 figure()
 hold on
 grid on
 plot(T,y)
 plot(T,-y)
 ylim([-.5 .5]);
 title('NACA 0012 Airfoil')
 xlabel('Chord Length [m]')
 ylabel('Thickness [m]')
 hold off


%% Functions

function [TotalArea] = CompTrapInt(panels,func)
        syms theta
        %func
        %class(func)
        %class(theta)
        panelWidth= (2*pi)/(panels);
        
        k = 0;
        k1 = k + panelWidth;
        
      %  fprintf('Trapezoidal Integration\n')
        
    for i = 1:panels
      %  fprintf('I am on interation: %i \n', i)
        func1 = subs(func,theta,k);
        func2 = subs(func,theta,k1);
        area(i) = panelWidth*(func1+func2)/2;
        k = k1;
        k1 = k+panelWidth;
        %area(i)
    end

    TotalArea = area;
end

function [Cl,Cd,Panels] = trap_Problem2(c,alpha,t,Cp_upper,Cp_lower)
% Inputs:  c = chord length, m
%          alpha = angle of attack, degrees
%          t = max airfoil thickness as a fraction of chord 
%          Cp_upper = Cp given for upper wall of airfoil, struct
%          Cp_lower = Cp given for lower wall of airfoil, struct
% Outputs: c_l = sectional coeff of lift (vector)
%          c_d = sectional coeff of drag (vector)
%          Panels = Number of panels used for each corresponding set of coeff. 
%                (vector)
% Methodology: To integrate the product of Cp and airfoil shape formula to
% solve for lift and drag coefficients per unit span

% Create a panel size vector to compare results of various fidelity
Panels = [1:5000 10000];
Panels = Panels';

% Set start/end for integration
a = 0;
b = c;

% Preallocate space for vectors
Cl = zeros(length(Panels),1);
Cd = zeros(length(Panels),1);
Cn = zeros(length(Panels),1);
Ca = zeros(length(Panels),1);

% Create for loop to calculate c_l and c_d for every panel size
disp('Starting for loop')

    for i = 1:length(Panels)
       N = Panels(i);

       % Create a x vector for every point along length of airfoil
       x = linspace(a,b,N+1);

       % Find C_p upper and lower for every x
       C_p_u = fnval(Cp_upper,x/c);
       C_p_l = fnval(Cp_lower,x/c);

       % Create dx for integration size
       panelWidth = x(2) - x(1);

       % Use body shape formula in integration
       y_t = (t/(0.2))*c*(0.2969*(sqrt(x./c)) - .1260*(x./c) - .3516*(x./c).^2 + .2843*(x./c).^3 - .1036*(x./c).^4);

       % Find dy_t/dx for integration product use
       for j = 2:length(y_t)
          dydx(j) = (y_t(j) - y_t(j-1))/panelWidth;  
       end

       % Create product to solve for force coefficients which will be integrated
       Fn = C_p_l - C_p_u;
       Fa = C_p_u.*dydx - C_p_l.*(-dydx);

       % Calculate c_n and c_a using trapezoidal rule with proper equations in
       % the integral, ignoring shear
       Cn(i,1) = (1/c)*(sum(Fn) - (Fn(1)) + (Fn(end)/2))*panelWidth;
       Ca(i,1) = (1/c)*(sum(Fa) - (Fa(1)) + (Fa(end)/2))*panelWidth; 

       % Calculate lift and drage coefficients from force coefficients and AoA
       Cl = Cn*cosd(alpha) - Ca*sind(alpha);
       Cd = Cn*sind(alpha) + Ca*cosd(alpha);
    end
end



function [TotalArea] = CompSimpInt(panels,func)
        syms theta
        %func
        %class(func)
        %class(theta)
        panelWidth= (2*pi)/(panels);
     
        
        k = 0;
        midpoint = k +panelWidth/2;
        k1 = k + panelWidth;
        
      %  fprintf('Simpson Integration\n')
        
    for i = 1:panels
       % fprintf('I am on interation: %i \n', i)
        func1 = subs(func,theta,k);
        func1_2 = subs(func,theta,midpoint);
        func2 = subs(func,theta,k1);
        area(i) = (panelWidth/6)*(func1+4*func1_2+func2);
        k = k1;
        midpoint = k+panelWidth/2;
        k1 = k+panelWidth;
        %area(i)
    end

    TotalArea = area;
end


