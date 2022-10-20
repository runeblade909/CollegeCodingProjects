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

% plotnum of which dataset to plot 1-18. If 0 no plots will be generated
plotnum = 18;

% True to plot Isp_SEM vs N Samples
plot_Isp = true;

%% Load Static Data

% File info
sample_freq = 1652; %[Hz]
pressure = 40; %[psi]
temp = 9; %[deg C]
m_water = 1; %[kg] Mass of water

% Find all files in folder, count them, then preallocate cell array
datafiles = dir('static_test_data\');
num_files = 18;
data = cell(1,num_files);

% Preallocate a num_files vec for use in many later preallocations
pre18 = zeros(1, num_files);

for i = 1:num_files
    % Add 3rd col of file data to second col of each cell
    tempdata = readmatrix(strcat('static_test_data\', datafiles(i+2).name));
    data{1,i}(:,2) = tempdata(:,3) * 4.44822; %[N]
    
    % Generate sample times (first col of each cell)
    data{1,i}(:,1) = (0:length(tempdata)-1)/sample_freq;

    % Each cell of 'data' looks like: (time [sec], thrust [N])
end

%% Integration Limits
% x1 & x2 are 1x18 vectors of times [in sec] of the start and end limits
% for integration and plotting 

x1 = [1.83 1.52 1.287 1.875 1.71 1.65 0.52 0.93 1.6 1.7 1.35 1.81 1.71 1.33 2.71 1.93 2.088 1.63];


peak = pre18; %[N]
peak_idx = pre18;  %[index]
peak_time = pre18; %[s] Time of peak from beg. of respective dataset
peak_time_adj = pre18; %[s] Time of peak from respective x1
slope = pre18;
x2 = pre18;
Last500 = pre18;
for s = 1:num_files
    % Peak thrusts & times
    [peak(1,s), peak_idx(1,s)] = max(data{1,s}(:,2));
    peak_time(1,s) = data{1,s}(peak_idx(1,s), 1);
    peak_time_adj(1,s) = peak_time(1,s) - x1(s);
   
end

for i = 1:num_files
Last500(1,i) = mean(data{1,i}(end-500:end,2));
    for k = peak_idx(1,i):length(data{1,i}(:,2))-51

                       slope(1,i) = (data{1,i}(k+150,2) - data {1,i}(k,2)) / (data{1,i}(k+150,1) - data {1,i}(k,1));
                       
                       if (slope(1,i) >= -1 && slope(1,i) <= 1) 
                           
                           x2(i) = data{1,i}(k+100,1); 
                           break
                       end
                       
    end
    
end
    


% x2 = x1 + 0.35;

%% Plot specified dataset and integration limits

% Plot the data specified by plotnum, x1, x2
if ~(plotnum == 0)
    figure
    hold on
    
    plot(data{1,plotnum}(:,1), data{1,plotnum}(:,2));
    yline(0, 'r', 'LineWidth', 0.75);
    xline(x1(plotnum), 'r', 'LineWidth', 0.75);
    xline(x2(plotnum), 'r', 'LineWidth', 0.75);
    
    xlim([x1(plotnum) - 0.5, x2(plotnum) + 0.5]);
    ylim([-2, max(data{1,plotnum}(:,2))+5]);
    
    title(strrep(datafiles(plotnum+2).name, '_', ' '));
    xlabel('Time [sec]');
    ylabel('Thrust [N]');
    hold off
end


%% Impulse & Specific Impulse
% Preallocate vectors
force_domain = cell(1, num_files);
I = pre18;

for n = 1:num_files
    time_step = data{1,n}(2,1); %[s] Time between measurements
    
    % Index of x1 and x2
    x1_idx = find(data{1,n}(:,1) >= x1(n), 1);
    x2_idx = find(data{1,n}(:,1) >= x2(n), 1);
    
    % Isolate forces from x1 to x2 then integrate to find impulse
    force_domain{1,n} = data{1,n}((x1_idx:x2_idx),2); %[N]
    I(1,n) = time_step * trapz(force_domain{1,n}); %[Ns] Impulse

end

Isp = I / (m_water * g0); %[s] Specific Impulse

%% Statistics
% Total time of thrust used for I
thrust_time = x2 - x1;

% Preallocate Vectors
% peak = pre18; %[N]
% peak_idx = pre18;  %[index]
% peak_time = pre18; %[s] Time of peak from beg. of respective dataset
% peak_time_adj = pre18; %[s] Time of peak from respective x1
Isp_mean = pre18; % Index is the average Isp for the first idx samples
Isp_std = pre18; %
Isp_SEM = pre18; % Index is the Isp SEM for first idx samples
Isp_confidence = zeros(3, num_files); % +/- Confidence ints for Isp (95%; 97.5%; 99%)

for s = 1:num_files
%     % Peak thrusts & times
%     [peak(1,s), peak_idx(1,s)] = max(data{1,s}(:,2));
%     peak_time(1,s) = data{1,s}(peak_idx(1,s), 1);
%     peak_time_adj(1,s) = peak_time(1,s) - x1(s);
    
    % Isp Std. Dev. & SEM
    Isp_mean(1,s) = mean(Isp(1,1:s));
    Isp_std(1,s) = std(Isp(1, 1:s));
    Isp_SEM(1,s) = Isp_std(1,s) / sqrt(s);
    Isp_confidence(:, s) = [1.96 * Isp_SEM(1,s); 2.24 * Isp_SEM(1,s); 2.58 * Isp_SEM(1,s)];
end

% Standard Deviation for peak thrusts and times
peak_mean = mean(peak); %[N] average peak thrust
peak_std = std(peak); %[N] Std. Dev. of 'peak'

peak_time_adj_mean = mean(peak_time_adj); %[s] Average peak time adj
peak_time_adj_std = std(peak_time_adj); %[s] Std. Dev. of 'peak_time_adj'

% Plot Isp_SEM vs N
if plot_Isp == true
    figure
    hold on
    plot(1:num_files, Isp_SEM);
    grid on
    title('Isp SEM vs N Samples');
    xlabel('N samples');
    ylabel('Isp SEM');
    hold off
end

