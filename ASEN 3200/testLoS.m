function inLoS = testLoS(r_site,r_sc,elevation_limit)
%DESCRIPTION: Determines whether the spacecraft is within line-of-sight
%(LoS) of the site given an elevation limit
%
%INPUT:
% r_site            The position vector of the site (km, 3x1)
% r_sc              The position vector of the spacecraft (km, 3x1)
% elevation_limit   Lower elevation limit (above the horizon) (rad)
%
%OUTPUT: 
% inLoS             A boolean flag (0 or 1); 1 indicates the spacecraft and
%                   the site have line-of-sight



%1) Compute whether the site and spacecraft have line of sight (hint, I
%suggest drawing a picture and writing this constraint as an inequality
%using a dot product)

r_sc_minus_site = r_sc - r_site; 
r_hat_sc_minus_site = r_sc_minus_site/norm(r_sc_minus_site); 
r_hat_site = r_site/norm(r_site); 
theta2 = acos(dot(r_hat_sc_minus_site, r_hat_site)); 
theta_horizon = pi/2 - theta2; 
inLoS = double(theta_horizon>=elevation_limit); 
end