%json file reading function: 
function [num_launches, num_spacecraft, satellite_list] = loadConstellation(fileName)
%DESCRIPTOIN: Ingests constellation description .json file and parses it
%into a list of structs with full initial orbit elements (km, s, rad) and
%satellite name.
%
%INPUTS:
% filename      A string indicating the name of the .json file to be parsed
%
%OUTPUTS:
% nl            Number of total launches
% ns            Total number of spacecraft between all launches
% satlist       Array of structs with 'name' and 'oe0' properties


%Temporary - just so the function runs the first time you use it.
%You'll need to change all of these!
num_launches = 0;
num_spacecraft = 0;
satellite_list.name = '';
satellite_list.oe0 = NaN(6,1);

%1) extract the constellation structure from the json file
openFile = fopen(fileName); 
rawData = fread(openFile,inf); 
charData = char(rawData'); 
fclose(openFile); 
data = jsondecode(charData); 
%2) read all of the launches and payloads to understand how many launches
% and spacecraft are in the constellation; note, this will be useful in
% Part 2!
num_launches = length(data.launches);

for i = 1:num_launches
    num_spacecraft = num_spacecraft + length(data.launches(i).payload);
end

    
%3) RECOMMENDED: Pre-allocate the satellite_list struct
satellite_list(num_spacecraft).name = '';
satellite_list(num_spacecraft).oe0 = NaN(6,1);

%4) Populate each entry in the satellite struct list with its name and
%initial orbit elements [a,e,i,Om,om,f] at time t0
count = 0; 
for i = 1:num_launches
    num_pay_per_launch = length(data.launches(i).payload);
    oe = [data.launches(i).orbit.a; data.launches(i).orbit.e; data.launches(i).orbit.i; data.launches(i).orbit.Om; data.launches(i).orbit.om; NaN]; 
    for j = 1:num_pay_per_launch
        count = count+1;
        satellite_list(count).name = data.launches(i).payload(j).name;
        f = data.launches(i).payload(j).f; 
        oe(end) = f; 
        satellite_list(count).oe0 = oe; 
    end
end

end

