%{

static_test.m
Bennett Grow, Zachary Reichenbach, Spring 2021, CU Boulder ASEN 2004 Bottle Rocket Lab

Loads static test stand data from all files in folder 'static_test_data'
into a cell array 'data' then calculates specific impulse for each dataset 
and various statistics

%}

clc
clear
close all
g0 = 9.81; %[m/s^2]




%% Load Static Data

% File info
sample_freq = 1652; %[Hz]
pressure = 40; %[psi]
temp = 11; %[deg C]
m_water = .6; %[kg] Mass of water

% Find all files in folder, count them, then preallocate cell array
datafiles = readtable('LA_Demo_S013_test2_600g');

data1 = table2array(datafiles(4:end,3));
    data = [((0:length(data1)-1)./sample_freq)' data1] ;

    data(4963,2) = (data(4963,2)*-1); 
    data(:,2) = data(:,2)*4.448;
%% Integration Limits
% x1 & x2 are 1x18 vectors of times [in sec] of the start and end limits
% for integration and plotting 

for k = 1:length(data(:,2))

                       slope = (data(k+100,2) - data (k,2)) / (data(k+100,1) - data(k,1));
                       
                       if (slope >= 100) 
                           x1 = data(k+100,1); 
                           break
                       end
                       
    end

% Preallocate a num_files vec for use in many later preallocations
pre18 = zeros(1, 1);
peak = pre18; %[N]
peak_idx = pre18;  %[index]
peak_time = pre18; %[s] Time of peak from beg. of respective dataset
peak_time_adj = pre18; %[s] Time of peak from respective x1
slope = pre18;
x2 = pre18;
Last500 = pre18;

    % Peak thrusts & times
    [peak, peak_idx] = max(data(:,2));
    peak_time = data(peak_idx, 1);
    peak_time_adj = peak_time - x1;



Last500 = mean(data(end-500:end,2));
    for k = peak_idx:length(data(:,2))-51

                       slope = (data(k+250,2) - data (k,2)) / (data(k+250,1) - data(k,1));
                       
                       if (slope >= -1 && slope <= 1) 
                           
                           x2 = data(k+100,1); 
                           break
                       end
                       
    end

    
lb= 4.448; %N

% x2 = x1 + 0.35;

%% Plot specified dataset and integration limits

% Plot the data specified by plotnum, x1, x2
    figure
    hold on
    
    plot(data(:,1), data(:,2));
    yline(0, 'r', 'LineWidth', 0.75);
    xline(x1, 'r', 'LineWidth', 0.75);
    xline(x2, 'r', 'LineWidth', 0.75);
    
    xlim([x1 - 0.5, x2 + 0.5]);
    ylim([-2, max(data(:,2))+5]);
    
    title('Isp of a 600g Propellant Bottle Rocket')
    xlabel('Time [sec]');
    ylabel('Thrust [N]');
    hold off



%% Impulse & Specific Impulse
% Preallocate vectors
force_domain = cell(1,1);
I = pre18;


    time_step = data(2,1); %[s] Time between measurements
    
    % Index of x1 and x2
    x1_idx = find(data(:,1) >= x1, 1);
    x2_idx = find(data(:,1) >= x2, 1);
    
    % Isolate forces from x1 to x2 then integrate to find impulse
    force_domain = data((x1_idx:x2_idx),2); %[N]
    I = time_step * trapz(force_domain); %[Ns] Impulse



Isp = I / (m_water * g0); %[s] Specific Impulse