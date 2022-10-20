function plot12(t,s)
    f = figure;
    f.Position = [50,50,900,700];
    subplot(4,3,1), plot(t,s(:,1));
    title('X Position (North +)');
    xlabel('Time [seconds]');
    ylabel('X Position [m]');
    subplot(4,3,2), plot(t,s(:,2));
    title('Y Position (East +)');
    xlabel('Time [seconds]');
    ylabel('Y Position [m]');
    subplot(4,3,3), plot(t,s(:,3));
    title('Z Position (Down +)');
    xlabel('Time [seconds]');
    ylabel('Z Position [m]');
    
    subplot(4,3,4), plot(t,s(:,4));
    title('Roll(\phi)');
    xlabel('Time [seconds]');
    ylabel('Roll(\phi) [rad]');
    subplot(4,3,5), plot(t,s(:,5));
    title('Pitch(\theta)');
    xlabel('Time [seconds]');
    ylabel('Pitch(\theta) [m]');
    subplot(4,3,6), plot(t,s(:,6));
    title('Yaw(\psi)');
    xlabel('Time [seconds]');
    ylabel('Yaw(\psi)Position [m]');
    
    subplot(4,3,7), plot(t,s(:,7));
    title('U Inertial Velocity (North +)');
    xlabel('Time [seconds]');
    ylabel('U Inertial Velocity [m]');
    subplot(4,3,8), plot(t,s(:,8));
    title('V Velocity (East +)');
    xlabel('Time [seconds]');
    ylabel('V Inertial Velocity [m]');
    subplot(4,3,9), plot(t,s(:,9));
    title('W Velocity (Down +)');
    xlabel('Time [seconds]');
    ylabel('W Inertial Velocity [m]');
    
    subplot(4,3,10), plot(t,s(:,10));
    title('p Angular Velocity');
    xlabel('Time [seconds]');
    ylabel('q Angular Velocity [m]');
    subplot(4,3,11), plot(t,s(:,11));
    title('p Angular Velocity');
    xlabel('Time [seconds]');
    ylabel('p Angular Velocity [m]');
    subplot(4,3,12), plot(t,s(:,12));
    title('r Angular Velocity');
    xlabel('Time [seconds]');
    ylabel('r Angular Velocity [m]');
    
    sgtitle('State Variables plots');
end