function x = odefun_level(t,s)
   eq_1 = firstEq(s(4),s(5),s(6),s(7),s(8),s(9));
   eq_2 = secondEq(s(4),s(5),s(6),s(10),s(11),s(12));
   eq_3 = [s(7);s(8);s(9)];
   eq_4 = [s(10);s(11);s(12)]; 
   
   x = [eq_1;eq_2;eq_3;eq_4];
end
function out = firstEq(phi,theta,psi,uE,vE,wE)
    A = [cosd(theta)*cosd(psi) sind(phi)*sind(theta)*cosd(psi)-cosd(phi)*sind(psi) cosd(phi)*sind(theta)*cosd(psi)+sind(phi)*sind(psi);...
        cosd(theta)*sind(psi) sind(phi)*sind(theta)*sind(psi)+cosd(phi)*cosd(psi) cosd(phi)*sind(theta)*sind(psi)-sind(phi)*cosd(psi);...
        -sind(theta) sind(phi)*cosd(theta) cosd(phi)*cosd(theta)];
    out = A*[uE;vE;wE];
end
function out = secondEq(phi,theta,psi,p,q,r)
    A = [1 sind(phi)*tand(theta) cosd(phi)*tand(theta);...
        0 cosd(phi) -sind(phi);...
        0 sind(phi)*secd(theta) cosd(phi)*secd(theta)];
    out = A*[p;q;r];
end
