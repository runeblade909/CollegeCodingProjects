clear; clc; close all;

% constellation = struct([]);
% launches = struct([]);
% orbit = struct([]);
% payload = struct([]);
constellation = struct();
launches = struct('orbit',{},'payload',{},'launchName',{});
constellation.companyName = 'Ethical Voluntary Idyllic Launcher';
constellation.companyExecutives = cell(3,1);
constellation.companyExecutives{1} = "Chris";
constellation.companyExecutives{2} = "Zak";
constellation.companyExecutives{3} = "Jon";
Numlaunches = 7; 
satellites = 350;
Re = 6378; %[km]
inc = pi/180 * - 50; %[rad]
e = .00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001; %Cirle
a = 1100+Re; %[km]
Omega = 0;
omega = 0;
k = 1;
Azimuth = 125;


for i = 1:Numlaunches
    f = 0;
    orbit.a = a;
    orbit.e = e;
    orbit.i = inc;
    orbit.Om = Omega;
    orbit.om = omega;
    orbit.AzimuthShift = Azimuth;
    for j = 1:satellites/Numlaunches
        payload(j).f = f;
        payload(j).name = num2str(k);
        f = f + 2*pi/(satellites/Numlaunches);
        k = k+1;
    end
    launches(i).orbit = orbit;
    launches(i).payload = payload;
    launches(i).launchName = num2str(i);
    Azimuth = Azimuth + 2*pi/Numlaunches;
    a = a-.1;
end
constellation.launches = launches;
output = jsonencode(constellation);
fid = fopen('Constellation.json','w');
fprintf(fid,output);
fclose(fid);