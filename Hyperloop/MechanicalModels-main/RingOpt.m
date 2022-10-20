%%
clc; clear; close all;

%%
rings_arr = 10:80;
arrsz = length(rings_arr);
rod_diam = 0.0254*(3/8); % 3/8" rod diameter
tarp_1ring_folded = 0.0254/2; % One section of tarp folded

for i = 1:arrsz

    rings = rings_arr(i);
    
    spacing = 31.9/(rings-1);

    len_compr = (rod_diam+tarp_1ring_folded)*(rings-2);

    len_arr(i) = len_compr+spacing;
end

figure(1);
plot(rings_arr, len_arr);
title("Length of back section (whole tunnel compressed + 1 ext section")
xlabel("Number of rings")
ylabel("Length of back section [m]")

%%
d_ring = 0.59;
h_max = d_ring-0.5;
t = 0.0254*0.01;
E = 324.5e6; % Young's mod from tensile test

w = 0.5; % 0.5 m width of tarp section
P = 28000;
Surf_A = w*t;
sensor = 9e-3;

sig_yield = 13463900; % Pa, from tensile test

% Pre-allocating all arrays to same size
a_mat = zeros(1,arrsz);
def = a_mat; strain = a_mat; Stress = a_mat; Stress2 = a_mat;

for i = 1:arrsz
    
    rings = rings_arr(i);
    spc = 31.9/(rings-1);
    
    syms a
    a = double(solve(h_max == a*cosh(spc/(2*a)) - a));
    a_mat(i) = a;
    
    syms x;
    diff_CAT = a*sinh(x/a);
    L_new = double(int(sqrt(1+diff_CAT^2), -spc/2, spc/2));

    def(i) = (L_new-spc)+sensor;
    strain(i) = def(i)/spc;
    
    Stress(i) = E*strain(i);

    Tension = P*t*a;
    Stress2(i) = Tension/Surf_A;

end

figure;
plot(rings_arr,Stress); hold on;
plot([rings_arr(1) rings_arr(end)], [sig_yield sig_yield],'r')
title('Stress in tarp - worst case scenario');
xlabel('Number of rings');
ylabel('Stress in tarp [Pa]');
legend('Main data', 'Yield stress of tarp')

figure;
Force = Stress * t * pi * d_ring;
plot(rings_arr, Force);
title('Force on tarp - worst case scenario');
xlabel('Number of rings');
ylabel('Force applied on tunnel [N]');

figure;
plot(rings_arr, Stress2, 'LineWidth', 2);
title('Actual stress in tarp based on catenary model');
xlabel('Number of rings');
ylabel('Stress in tarp [Pa]');

%% 

figure(4); hold on;

N = 25;
colorss = ['k','b','r','m'];

for i = 0:3
    index = i*20 + 1;
    
    a = a_mat(index);
    
    spc = 31.9/(rings_arr(index)-1);
    x_arr = linspace(-spc/2, spc/2, N);
    y_arr = a.*cosh(x_arr);
    
    maxim = max(y_arr);
    y_arr = (y_arr-maxim) + 0.59;

    plot(x_arr, y_arr, colorss(i+1));
end


% 45 rings
index = 37;
a = a_mat(index);
spc = 31.9/(rings_arr(index)-1);
x_arr = linspace(-spc/2, spc/2, N);
y_arr = a.*cosh(x_arr);
maxim = max(y_arr);
y_arr = (y_arr-maxim) + 0.59;
plot(x_arr, y_arr, 'linewidth', 2);

plot([-2 2],[0.5 0.5],'--r');

ylim([0 0.6]);
legend('10 rings', '30 rings', '50 rings', '70 rings', '**46 rings**');


xlabel("Length along one section [meters]");
ylabel("Tarp sag [meters]");
title("Visualizing tarp sag");

%% FOS

figure;
FOS = sig_yield./Stress2;
plot(rings_arr, FOS);
