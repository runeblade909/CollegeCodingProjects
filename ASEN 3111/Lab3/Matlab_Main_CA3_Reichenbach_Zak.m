% Computational Lab #4: Prandtl Lifting Line Theory
%Zak Reichenbach
%10/25/2021

%% Housekeeping

close all; clear all; clc;

%% Problem 1: Implement Prandtl Lifting Line Theory Code

%First Verification Check

b = 80; 
a0_t = 2*pi; 
a0_r = 2*pi; 
c_t = 10; 
c_r = 10; 
aero_t = 0; 
aero_r = 0; 
geo_t = 5; 
geo_r = 5; 
N = 6;

[e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);

fprintf('PROBLEM 1 FIRST VERIFICATION TEST\n');
fprintf('  e: %f \n',e);
fprintf('  c_L: %f \n',c_L);
fprintf('  c_Di: %f \n',c_Di);
fprintf('\n');

%Second Verification Check (Slightly tapered wings)

b = 80; 
a0_t = 2*pi; 
a0_r = 2*pi; 
c_t = 8; 
c_r = 10; 
aero_t = 0; 
aero_r = 0; 
geo_t = 5; 
geo_r = 5; 
N = 6;
[e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);

fprintf('PROBLEM 1 FIRST VERIFICATION TEST 2\n');
fprintf('  e: %f \n',e);
fprintf('  c_L: %f \n',c_L);
fprintf('  c_Di: %f \n',c_Di);
fprintf('\n');

%Third Verification Check (Geometric Twist)

b = 80; 
a0_t = 2*pi; 
a0_r = 2*pi; 
c_t = 10; 
c_r = 10; 
aero_t = 0; 
aero_r = 0; 
geo_t = 0; 
geo_r = 5; 
N = 6;
[e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);

fprintf('PROBLEM 1 FIRST VERIFICATION TEST 3\n');
fprintf('  e: %f \n',e);
fprintf('  c_L: %f \n',c_L);
fprintf('  c_Di: %f \n',c_Di);
fprintf('\n');


%Fourth Verification Check (Aspec Ratio B^2/A)

b = 100; 
a0_t = 2*pi; 
a0_r = 2*pi; 
c_t = 10; 
c_r = 10; 
aero_t = 0; 
aero_r = 0; 
geo_t = 5; 
geo_r = 5; 
N = 6;
[e,c_L,c_Di] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);

fprintf('PROBLEM 1 FIRST VERIFICATION TEST 4\n');
fprintf('  e: %f \n',e);
fprintf('  c_L: %f \n',c_L);
fprintf('  c_Di: %f \n',c_Di);
fprintf('\n');




%% Problem 2: Apply Prandtl Lifting Line Theory Code

b = 80;
c_t = 4;
c_r = 12;
geo_t = 1;
geo_r = 6;

V_inf = 130*5280/3600; %ft/sec
rho_inf = 0.002377; %slug/ft
S = (c_t+c_r)*b/2;  %Span

%Airfoil Dimensions

[XB_0012, YB_0012] = NACA_Airfoil(0,0,0.12*c_t,c_t,100);
[XB_2412, YB_2412] = NACA_Airfoil(0.02,0.4,0.12*c_r,c_r,100);

alpha = linspace(-5,10,10);


%Vortex Panel Method
for i = 1:length(alpha)
    [cl_0012(i),~,~] = Vortex_Panel(XB_0012,YB_0012,V_inf,alpha(i),0);
    [cl_2412(i),~,~] = Vortex_Panel(XB_2412,YB_2412,V_inf,alpha(i),0);
end

%

lfit_0012 = polyfit(alpha,cl_0012,1);
a0_t = lfit_0012(1)*180/pi;
aero_t = (-lfit_0012(2)./lfit_0012(1));

lfit_2412 = polyfit(alpha,cl_2412,1);
a0_r = lfit_2412(1)*180/pi;
aero_r = (-lfit_2412(2)./lfit_2412(1));

%PLLT Calculations

N = 300;

[e_exact,c_L_exact,c_Di_exact] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);

L_exact = 0.5*rho_inf*V_inf^2*S*c_L_exact;
Di_exact = 0.5*rho_inf*V_inf^2*S*c_Di_exact;

fprintf('PROBLEM 2 EXACT ANSWERS (USING 300 ODD TERMS)\n');
fprintf('  e: %f \n',e_exact);
fprintf('  c_L: %f \n',c_L_exact);
fprintf('  c_Di: %f \n',c_Di_exact);
fprintf('  L: %f pounds\n',L_exact);
fprintf('  Di: %f pounds\n',Di_exact);
fprintf('\n');

for N = 1:300
    [e(N), c_L(N), c_Di(N)] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);
end

%Find Relative Error

N = 1:300;

E_L = 100*abs(c_L - c_L_exact)/c_L_exact;
E_Di = 100*abs(c_Di - c_Di_exact)/c_Di_exact;

%5 Percent Relative Error

EL5 = lt(E_L,5);
i5 = find(EL5, 1, 'first');
NL5 = N(i5);
EDi5 = lt(E_Di,5);
i5 = find(EDi5, 1, 'first');
NDi5 = N(i5);

%2 Percent Relative Error

EL2 = lt(E_L,2);
i2 = find(EL2, 1, 'first');
NL2 = N(i2);
EDi2 = lt(E_Di,2);
i2 = find(EDi2, 1, 'first');
NDi1 = N(i2);

% 2/10 Percent Relative Error

EL02 = lt(E_L, 0.2);
i02 = find(EL02, 1, 'first');
NL02 = N(i02);
EDi02 = lt(E_Di, 0.2);
i02 = find(EDi02, 1, 'first');
NDi01 = N(i02);

%Window Response

fprintf('PROBLEM 2 LIFT CONVERGENCE \n');
fprintf('  The lift converges to within 5 percent using %i panels.\n', NL5);
fprintf('  The lift converges to within 2 percent using %i panels.\n', NL2);
fprintf('  The lift converges to within 0.2 percent using %i panels.\n', NL02);
fprintf('\n');

fprintf('PROBLEM 2 INDUCED DRAG CONVERGENCE \n');
fprintf('  The induced drag converges to within 5 percent using %i panels.\n', NDi5);
fprintf('  The induced drag converges to within 2 percent using %i panels.\n', NDi1);
fprintf('  The induced drag converges to within 0.2 percent using %i panels.\n', NDi01);
fprintf('\n');

%% Problem 3: Study on Taper Ratio/Aspect Ratio Effect on Span Effiency Factor

%Wing Properties
c_r = 15; 
a0_t = 2*pi; 
a0_r = 2*pi; 
aero_t = 0; 
aero_r = 0; 
geo_t = 5; 
geo_r = 5; 
N = 100;

TR = 0:0.02:1;
AR = [4 6 8 10];

%Varying Taper Ratio

for i = 1:size(TR,2)
    for j = 1:size(AR,2)
        c_t = TR(i)*c_r;
        b = 0.5*AR(j)*(c_t+c_r);
        
        [e(i,j), ~, ~] = PLLT(b,a0_t,a0_r,c_t,c_r,aero_t,aero_r,geo_t,geo_r,N);
    end
end

    figure()
    plot(TR,e(:,1),TR,e(:,2),TR,e(:,3),TR,e(:,4));
    set(gca,'FontName','Times New Roman','FontSize',14);
    title('Span Effiency Factor Versus Taper Ratio','FontName','Times New Roman','FontSize',16);
    xlabel('Taper Ratio $c_t/c_r$','Interpreter','Latex');
    ylabel('Span Effiency Factor $e$','Interpreter','Latex');
    l = legend('$AR = 4$','$AR = 6$','$AR = 8$','$AR = 10$');
    set(l,'Interpreter','Latex','Location','Best','FontSize',16);
    
    