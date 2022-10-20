%orbit state derivative for ode45 simulation:
function dXdt = orbit(mu,X)
r = X(1:3);
rdot = X(4:6);
rddot = -mu/(norm(r))^3 * r; 
dXdt = [rdot; rddot];
end