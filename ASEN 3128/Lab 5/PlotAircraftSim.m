function PlotAircraftSim(TOUT, aircraft_state, control_surfaces, col, wind,part)


%%%%%%%%%%%%%%%%%%%%%%%%
%figure;
figure(1+part)
 
subplot(311);
h1= plot(TOUT, aircraft_state(:,1),col);hold on;
title('Position v Time');   
ylabel('X [m]')    

subplot(312);
plot(TOUT, aircraft_state(:,2),col);hold on;
 ylabel('Y [m]')    
 
subplot(313);
plot(TOUT, aircraft_state(:,3),col);hold on;
ylabel('Z [m]')    
xlabel('time [sec]');

%%%%%%%%%%%%%%%%%%%%%%%%
%figure;
figure(2+part)
 
subplot(311);
plot(TOUT, (180/pi)*aircraft_state(:,4),col);hold on;
title('Euler Angles v Time');   
ylabel('Roll [deg]')    

subplot(312);
plot(TOUT, (180/pi)*aircraft_state(:,5),col);hold on;
 ylabel('Pitch [deg]')    
 
subplot(313);
plot(TOUT, (180/pi)*aircraft_state(:,6),col);hold on;
ylabel('Yaw [deg]')    
xlabel('time [sec]');

%%%%%%%%%%%%%%%%%%%%%%%%
%figure;
figure(3+part)
 
subplot(311);
plot(TOUT, aircraft_state(:,7),col);hold on;
title('Velocity v Time');   
ylabel('uE [m/s]')    

subplot(312);
plot(TOUT, aircraft_state(:,8),col);hold on;
 ylabel('vE [m/s]')    
 
subplot(313);
plot(TOUT, aircraft_state(:,9),col);hold on;
ylabel('wE [m/s]')    
xlabel('time [sec]');

%%%%%%%%%%%%%%%%%%%%%%%%
%figure;
figure(4+part)
 
subplot(311);
plot(TOUT, (180/pi)*aircraft_state(:,10),col);hold on;
title('Angular Velocity v Time');   
ylabel('p [deg/s]')    

subplot(312);
plot(TOUT, (180/pi)*aircraft_state(:,11),col);hold on;
 ylabel('q [deg/s]')    
 
subplot(313);
plot(TOUT, (180/pi)*aircraft_state(:,12),col);hold on;
ylabel('r [deg/s]')    
xlabel('time [sec]');

%%%%%%%%%%%%%%%%%%%%%%%%
%figure;
figure(5+part)
plot3(aircraft_state(:,1),aircraft_state(:,2),-aircraft_state(:,3),col);hold on;
title('Aircraft Flight Path'); 
ylabel('Y Position')    
xlabel('X Position');
zlabel('Z Position');

%%%%%%%%%%%%%%%%%%%%%%%%
if (~isempty(control_surfaces))
    %figure;
    figure(6+part)
    
    subplot(411);
    plot(TOUT, control_surfaces(:,1),col);hold on;
    title('Control Surfaces v Time');   
    ylabel('Elevator [rad]')    

    subplot(412);
    plot(TOUT, control_surfaces(:,2),col);hold on;
    ylabel('Aileron [rad]')      
 
    subplot(413);
    plot(TOUT, control_surfaces(:,3),col);hold on;
    ylabel('Rudder [rad]')       
    
    subplot(414);
    plot(TOUT, control_surfaces(:,4),col);hold on;
    ylabel('Throttle [frac]')     
    xlabel('time [sec]');
end

    %figure;
    figure(7+part)
    
    subplot(311);
    plot(TOUT, wind(:,1),col);hold on;
    title('Euler Angles for Wind v Time');   
    ylabel('V (m/s)')    