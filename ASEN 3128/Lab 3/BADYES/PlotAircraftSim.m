function PlotAircraftSim(t,state,control,fig,col,num)

    %% Plot position vs time
    figure(fig(1))
    subplot(3,1,1)
    plot(t,state(:,1),col,'linewidth',2); hold on;
    grid on
    ylabel('x');
    title([num ' Position vs. time']);
    subplot(3,1,2)
    plot(t,state(:,2),col,'linewidth',2); hold on;
    grid on
    ylabel('y');
    subplot(3,1,3)
    plot(t,state(:,3),col,'linewidth',2); hold on;
    grid on
    ylabel('z');
    xlabel('time [s]');
    
    %% Plot Euler andgles vs time
    figure(fig(2))
    subplot(3,1,1)
    plot(t,state(:,4),col,'linewidth',2); hold on;
    grid on
    ylabel('\phi');
    title([num ' Euler angles vs. time']);
    subplot(3,1,2)
    plot(t,state(:,5),col,'linewidth',2); hold on;
    grid on
    ylabel('\theta');
    subplot(3,1,3)
    plot(t,state(:,6),col,'linewidth',2); hold on;
    grid on
    ylabel('\psi');
    xlabel('time [s]');
    
    %% Plot velocity vs time
    figure(fig(3))
    subplot(3,1,1)
    plot(t,state(:,7),col,'linewidth',2); hold on;
    grid on
    ylabel('u');
    title([num ' Inertial Velocity vs. time']);
    subplot(3,1,2)
    plot(t,state(:,8),col,'linewidth',2); hold on;
    grid on
    ylabel('v');
    subplot(3,1,3)
    plot(t,state(:,9),col,'linewidth',2); hold on;
    grid on
    ylabel('w');
    xlabel('time [s]');
    
    %% Plot angular velocity vs time    
    figure(fig(4))
    subplot(3,1,1)
    plot(t,state(:,10),col,'linewidth',2); hold on;
    grid on
    ylabel('p');
    title([num ' Angular Velocity vs. time']);
    subplot(3,1,2)
    plot(t,state(:,11),col,'linewidth',2); hold on;
    grid on
    ylabel('q');
    subplot(3,1,3)
    plot(t,state(:,12),col,'linewidth',2); hold on;
    grid on
    ylabel('r');
    xlabel('time [s]');
    
    %% Plot control forces and moments vs time
    figure(fig(5))
    subplot(2,2,1)
    plot(t,control(:,1),col,'linewidth',2); hold on;
    grid on
    title([num ' Control Force vs. time']);
    ylabel('Control Force');
    xlim([t(1) t(end)]);
    subplot(2,2,2)
    plot(t,control(:,2),col,'linewidth',2); hold on;
    grid on
    title([num ' x Control Moment vs. time']);
    ylabel('x Control Moment');
    xlim([t(1) t(end)]);
    subplot(2,2,3)
    plot(t,control(:,3),col,'linewidth',2); hold on;
    grid on
    title([num ' y Control Moment vs. time']);
    ylabel('y Control Moment');
    xlabel('time [s]');
    xlim([t(1) t(end)]);
    subplot(2,2,4)
    plot(t,control(:,4),col,'linewidth',2); hold on;
    grid on
    title([num ' z Control Moment vs. time']);
    xlabel('time [s]');
    ylabel('z Control Moment');
    xlim([t(1) t(end)]);
    
    %% Plot flight path of the simulated quadrotor
    figure(fig(6))
    plot3(state(:,1),state(:,2),state(:,3),col,'linewidth',2); hold on;
    plot3(state(1,1),state(1,2),state(1,3),'g.','markersize',20); hold on;
    plot3(state(end,1),state(end,2),state(end,3),'r.','markersize',20); hold on;
    set(gca, 'YDir','reverse')
    set(gca, 'ZDir','reverse')
    grid on
    title([num ' Flight Path']);
    xlabel('x');
    ylabel('y');
    zlabel('z');

end