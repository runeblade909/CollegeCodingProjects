%% DataFetch
%
% ASEN 2002 Lab 2 Group 3
% Zak Reichenbach, Cate Billings, Amanda Marlow, & Bennett Grow
%
% Loads boundary, pitot, venturi, and manometer data into labdata.mat
% If the command window successfully clears the data is correctly loaded.
% 
% Examples of how to access a single cell of saved data:
%   boundarydata
%       boundarydata{1}.BoundaryLayer_S301_1{1,1}
%       where the first number is the port number
%   pitotdata
%       pitotdata.VelocityVoltage_S301_1{1,1}
%   venturidata
%       venturidata.VelocityVoltage_S301_3{1,1}
%   manometerdata
%       standard table
%       manometerdata{1,1}
%  
%% Import from 2019 > Boundary
%files_boundaryLayer = dir('2002 Lab 2/ASEN 2002/Aero_Lab_1_2019_Group_Data/BoundaryLayerData/');
files_boundaryLayer = dir('Aero_Lab_1_2019_Group_Data/BoundaryLayerData/');

clc;

outerrange = 4:14;
n = length(files_boundaryLayer);
namelist = cell(1, n);
valuelist = cell(1, n);
vlist = 1:n;
boundarydata = cell(1, 11);
nameboundary = cell(1, 11);
a = struct();


for i = 1:n
    if (i >= outerrange(1)) && (i <= max(outerrange))
        % Add portnames
        %namelist(i) = cellstr(files_boundaryLayer(i).name);
        name = files_boundaryLayer(i).name;
        %portnum = str2double(name(end));
        portnum = sscanf(name, 'Port_%d');
        
        % Allocate lists for the contents of a port folder
        files_port = dir(strcat(files_boundaryLayer(i).folder,'/',files_boundaryLayer(i).name));
        tempvaluelist = cell(1, length(files_port)-2);
        tempnamelist = cell(1, length(files_port)-2);
        
        % Loop through contents of a port folder
        for j = 3:length(files_port)
            % Record names of files
            [~, tempname, ~] = fileparts(files_port(j).name);
            tempnamelist(j-2) = cellstr(tempname);
            
            % Record data from files
            tempfile = load(strcat(files_port(j).folder,'/',files_port(j).name));
            tempvaluelist{j-2} = num2cell(tempfile);
        end
        nameboundary{portnum} = tempnamelist;
        
        % Put together names and data from files into a struct
        tempstruct = cell2struct(tempvaluelist,tempnamelist,2);
        
        % Add the struct to a cell array
        %boundarydata{portnum} = struct2cell(tempstruct);
        boundarydata{portnum} = tempstruct;
        % Access data like:
            % boundarydata{1}.BoundaryLayer_S301_1{1,1}
        
    end
end


%% Import from 2019 > Velocity > Pitot % Venturi
%files_VelocityVoltage_Pitot = dir('2002 Lab 2/ASEN 2002/Aero_Lab_1_2019_Group_Data/VelocityVoltageData/PitotProbeToPressureTransducer/');
%files_VelocityVoltage_Venturi = dir('2002 Lab 2/ASEN 2002/Aero_Lab_1_2019_Group_Data/VelocityVoltageData/VenturiTubeToPressureTransducer/');
files_VelocityVoltage_Pitot = dir('Aero_Lab_1_2019_Group_Data/VelocityVoltageData/PitotProbeToPressureTransducer/');
files_VelocityVoltage_Venturi = dir('Aero_Lab_1_2019_Group_Data/VelocityVoltageData/VenturiTubeToPressureTransducer/');
range = 3:14; % Specify the range of useful files
n = length(files_VelocityVoltage_Pitot); % number of 'files'
% Preallocate cellarrays for fields and values of struct that will store file data
fieldlist_Pitot = cell(1, n);
valuelist_Pitot = cell(1, n);

fieldlist_Venturi = cell(1, n);
valuelist_Venturi = cell(1, n);

namepitot = strings(1,12);
nameventuri = strings(1,12);
j = 1;

% Loop through files in 2019 > Velocity > Pitot/Venturi
for i = 1:n
% Store filenames without .csv for struct fields
    [~, tempnamepitot, ~] = fileparts(files_VelocityVoltage_Pitot(i).name);
    fieldlist_Pitot(i) = cellstr(tempnamepitot);
    
    
    [~, tempnameventuri, ~] = fileparts(files_VelocityVoltage_Venturi(i).name);
    fieldlist_Venturi(i) = cellstr(tempnameventuri);
    
    %valuelist(i) = cellstr(strcat("cell",int2str(i)));
    
    % Store data from each file in struct
    if (i >= range(1)) && (i <= max(range)) % Only access files from range
        % Load data to be stored from file
        tempfile = load(strcat(files_VelocityVoltage_Pitot(i).folder,'/',files_VelocityVoltage_Pitot(i).name));
        % Store the data in a cell
        valuelist_Pitot{i} = num2cell(tempfile);
        
        % Load data to be stored from file
        tempfile = load(strcat(files_VelocityVoltage_Venturi(i).folder,'/',files_VelocityVoltage_Venturi(i).name));
        % Store the data in a cell
        valuelist_Venturi{i} = num2cell(tempfile);
        
       %Names of the files
        namepitot(j) = tempnamepitot;
        nameventuri(j) = tempnameventuri;
        j=j+1;
    end
end
% Combine the fields and their values into a struct over the desired range
pitotdata = cell2struct(valuelist_Pitot(range), fieldlist_Pitot(range), 2);
venturidata = cell2struct(valuelist_Venturi(range), fieldlist_Venturi(range), 2);
%% Import manometer data
% files_VelocityVoltage_Manometer = dir('2002 Lab 2/ASEN 2002/Aero_Lab_1_2019_Group_Data/VelocityVoltageData/');
files_VelocityVoltage_Manometer = dir('Aero_Lab_1_2019_Group_Data/VelocityVoltageData/');
manometerdata = readtable(files_VelocityVoltage_Manometer(6).name);
%  = readmatrix(fullfile(pwd, "Aero_Lab_1_2019_Group_Data", "VelocityVoltageData", "manometer.xlsx"));

%% Import Airfoil Pressure Data
files_AirfoilPressure = dir('2002_Aero_Lab_2_Group_Data/');
airfoil_range = 3:28;
n = length(files_AirfoilPressure); % number of 'files'
% Preallocate cellarrays for fields and values of struct that will store file data
fieldlist_Airfoil = cell(1, n);
valuelist_Airfoil = cell(1, n);
nameairfoil = strings(1, length(airfoil_range));
j = 1;
for i = 1:n
    % Store filenames without .csv for struct fields
    [~, tempnameairfoil, ~] = fileparts(files_AirfoilPressure(i).name);
    fieldlist_Airfoil(i) = cellstr(tempnameairfoil);
    
     % Store data from each file in struct
    if (i >= airfoil_range(1)) && (i <= max(airfoil_range)) % Only access files from range
        % Load data to be stored from file
        tempfile = load(strcat(files_AirfoilPressure(i).folder,'/',files_AirfoilPressure(i).name));
        % Store the data in a cell
        valuelist_Airfoil{i} = num2cell(tempfile);
        
       %Names of the files
        nameairfoil(j) = tempnameairfoil;
        j=j+1;
    end
end
airfoildata = cell2struct(valuelist_Airfoil(airfoil_range), fieldlist_Airfoil(airfoil_range), 2);



%% Save workspace data
save labdata.mat boundarydata pitotdata venturidata manometerdata namepitot nameventuri nameboundary airfoildata nameairfoil;
clear
clc
load labdata.mat



