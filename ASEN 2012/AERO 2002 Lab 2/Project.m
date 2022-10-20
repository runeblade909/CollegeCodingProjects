% Group 3
%Zak Reichenbach
%Cate Billings45


%% House keeping
clear; 
clc;

%% Import from 2019 > Boundary

%files_boundaryLayer = dir('Aero_Lab_1_2019_Group_Data/BoundaryLayerData/');
%data_Boundary = strcat(files_boundaryLayer(4).folder,'/',files_boundaryLayer(4).name); 

%% Import from 2019 > Velocity > Pitot
files_VelocityVoltage_Pitot = dir('Aero Lab 1 - 2019 Group Data\VelocityVoltageData\PitotProbeToPressureTransducer');
range = 3:14; % Specify the range of useful files
n = length(files_VelocityVoltage_Pitot); % number of 'files'

% Preallocate cellarrays for fields and values of struct that will store file data
fieldlist = cell(1, n);
valuelist = cell(1, n);

% Loop through files in 2019 > Velocity > Pitot
for i = 1:n
    % Store filenames without .csv for struct fields
    [~, tempname, ~] = fileparts(files_VelocityVoltage_Pitot(i).name);
    fieldlist(i) = cellstr(tempname);
    
    %valuelist(i) = cellstr(strcat("cell",int2str(i)));
    
    % Store data from each file in struct
    if (i >= range(1)) && (i <= max(range)) % Only access files from range
        % Load data to be stored from file
        tempfile = load(strcat(files_VelocityVoltage_Pitot(i).folder,'/',files_VelocityVoltage_Pitot(i).name));
        % Store the data in a cell
        valuelist{i} = num2cell(tempfile);
    end
end

% Combine the fields and their values into a struct over the desired range
pitotdata = cell2struct(valuelist(range), fieldlist(range), 2);

%% Import from 2019 > Velocity > Venturi
% files_VelocityVoltage_Venturi = dir('Aero_Lab_1_2019_Group_Data/VelocityVoltageData/VenturiTubeToPressureTransducer/');
% data_Venturi = load(strcat(files_VelocityVoltage_Venturi(7).folder,'/',files_VelocityVoltage_Venturi(7).name)); %@302_3?

%% Import manometer data
% files_VelocityVoltage_Manometer = dir('Aero_Lab_1_2019_Group_Data/VelocityVoltageData/');
% data_Manometer = readtable(files_VelocityVoltage_Manometer(6).name);
manometerdata = readmatrix(fullfile(pwd, "Aero Lab 1 - 2019 Group Data", "VelocityVoltageData", "manometer.xlsx"));


%% Save workspace data
save('labdata')





%% OLD


% % Group 3
% % Zak Reichenbach
% % Cate Billings
% % Amanda Marlow
% % Bennett Grow
% 
% % House keeping
% clear; clc;
% 
% % Import Data ONLY
% 
% % files_boundaryLayer = dir('Aero_Lab_1_2019_Group_Data/BoundaryLayerData/');
% % sum = 0;
% % for i = 4:14
% %     data_Boundary = (strcat(files_boundaryLayer(i).folder,'/',files_boundaryLayer(i).name)); 
% %     for j = data_Boundary
% %         name = strcat('Port',i - 3,'File',j);
% %     end
% % end
% 
% 
% files_VelocityVoltage_Pitot = dir('Aero_Lab_1_2019_Group_Data/VelocityVoltageData/PitotProbeToPressureTransducer/');
% data_Pitot = load(strcat(files_VelocityVoltage_Pitot(4).folder,'/',files_VelocityVoltage_Pitot(4).name)); %@301_2?
% 
% files_VelocityVoltage_Venturi = dir('Aero_Lab_1_2019_Group_Data/VelocityVoltageData/VenturiTubeToPressureTransducer/');
% data_Venturi = load(strcat(files_VelocityVoltage_Venturi(7).folder,'/',files_VelocityVoltage_Venturi(7).name)); %@302_3?
% 
% files_VelocityVoltage_Manometer = dir('Aero_Lab_1_2019_Group_Data/VelocityVoltageData/');
% data_Manometer = readtable(files_VelocityVoltage_Manometer(6).name);% READ
% 
% % TABLE SEEMS TO WORK TO GET THE DATA FOR THE MANOMETER
% 
% 
% 
% manometerdata = readmatrix(fullfile(pwd, "Aero_Lab_1_2019_Group_Data", "VelocityVoltageData", "manometer.xlsx"));
% 
% 
% save('labdata')