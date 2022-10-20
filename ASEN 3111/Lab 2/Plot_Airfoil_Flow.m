function [x,y,Psi,Phi,P,levels,bounds,N_final,err_vec] = Plot_Airfoil_Flow(c,alpha,V_inf,P_inf,rho_inf,N,step)
% Function to calculate and plot the stream function, velocity potential,
% and the pressure contour for the given inputs
%
% Inputs:   c       - chord length
%           alpha   - angle of attack
%           V_inf   - free stream velocity
%           P_inf   - free stream pressure
%           rho_inf - free stream density
%           N       - number of vortecies
%
% Outputs:  Psi     - stream function
%           Phi     - velocity potential
%           P       - pressure contour
    % Author: Zak Reichenbach
    % Date: 10/11/2021


%% If step == 1
    if step == 1
%       1 - regular function that calculates Psi, Phi, and P for 50000
%           vortices and then finds the number of vortices to get less 
%           than a 0.5% error

        % Domain Declaration
        x_min = -c/2;
        x_max = c+c/2;
        y_min = -c/2;
        y_max = c/2;
        bounds = [x_min x_max y_min y_max];

        % Define Number of Grid Points
        n_x = 100; % steps in the x direction
        n_y = 100; % steps in the y direction

        % Define Mesh
        %Assume that 50000 is enough to be about 100% accurate
        N_2 = 50000;
        [x,y] = meshgrid(linspace(x_min,x_max,n_x),linspace(y_min,y_max,n_y));

        % Define radius
        r =@(x1) sqrt((x-x1).^2 + (y).^2);
        
        % Define Vorticies
        delta_x = c/N_2; %panel width
        x_vortex = linspace(delta_x/2,c-delta_x,N_2); %x-values along chord

        % Calculate vortex strength and circulation
        gamma = 2*alpha*V_inf*sqrt((1-(x_vortex/c))./(x_vortex/c)); % vortex strength
        Gamma = gamma*delta_x;                                      % circulation

        % Calculate psi for uniform flow (Eq. 3.55)
        Psi_uniform = V_inf*(y*cos(alpha)-x*sin(alpha));    

        % Calculate phi for uniform flow (Eq. 3.53)
        Phi_uniform = V_inf*(x*cos(alpha)-y*sin(alpha));    

        % Calculate psi for vortex (Eq. 3.114)
        Psi_vortex = 0;
        for i = 1:N_2
            % calculate the stream function of each vortex and add them
            % using super position to get the total stream function
            Psi_vortex = Psi_vortex + (Gamma(i)*log(r(x_vortex(i))))/(2*pi);
        end

        % Calulate phi for vortex (Eq. 3.112)
        Phi_vortex = 0;
        for i = 1:N_2
            % calculate the velocity potential of each vortices and add them
            % using super position to get the total velocity potential
            theta = atan2(-y,-x+x_vortex(i));
            Phi_vortex = Phi_vortex + -(Gamma(i)*theta)/(2*pi);
        end

        % Add all the stream functions
        Psi_N = Psi_vortex + Psi_uniform;

        % Add all the velocity potentials
        Phi_N = Phi_vortex + Phi_uniform; 

        % Calculate P
        q_inf = (1/2)*rho_inf*V_inf^2; % dynamic pressure
        Cp = 1-(gradient(Phi_N,(x_max-x_min)/n_x,(y_max-y_min)/n_y)./V_inf).^2; %  Cp = 1 - ([Vx Vy]/V)^2 from cylindrical flow
        % Add all the pressure
        P_N = P_inf + Cp*q_inf;

        % Initialize error
        err = 1; 
        j = 1;
        err_vec = zeros(length(N:100:24000),1);

        % While loop to find the number of vortices needed for error less than 0.5
        while err > 0.005

            % Define Mesh
            [x,y] = meshgrid(linspace(x_min,x_max,n_x),linspace(y_min,y_max,n_y));

            % Define radius
            r =@(x1) sqrt((x-x1).^2 + (y).^2);

            % Define vortcies
            delta_x = c/N;
            x_vortex = linspace(delta_x/2,c-delta_x,N);

            % Calculate vortex strength and circulation
            gamma = 2*alpha*V_inf*sqrt((1-(x_vortex/c))./(x_vortex/c)); % vortex strength
            Gamma = gamma*delta_x;                                      % circulation

            % Calculate psi for uniform stream (Eq. 3.55)
            Psi_uniform = V_inf*(y*cos(alpha)-x*sin(alpha));

            % Calculate phi for uniform stream (Eq. 3.53)
            Phi_uniform = V_inf*(x*cos(alpha)-y*sin(alpha));

            % Calculate psi for vortex (Eq. 3.114)
            Psi_vortex = 0;
            for i = 1:N
                % calculate the stream function of each vortices and add them
                % together to get the total stream function of the system
                Psi_vortex = Psi_vortex + (Gamma(i)*log(r(x_vortex(i))))/(2*pi);
            end

            % Calulate phi for vortex (Eq. 3.112)
            Phi_vortex = 0;
            for i = 1:N
                % calculate the velocity potential of each vortices and add them
                % together to get the total velocity potential of the system
                theta = atan2(-y,-x+x_vortex(i));
                Phi_vortex = Phi_vortex + -(Gamma(i)*theta)/(2*pi);
            end

            % Add all the stream functions
            Psi = Psi_vortex + Psi_uniform;

            % Add all the velocity potentials
            Phi = Phi_vortex + Phi_uniform;

            % Calculate P
            q_inf = (1/2)*rho_inf*V_inf^2; % dynamic pressure
            Cp = 1-(gradient(Phi,(x_max-x_min)/n_x,(y_max-y_min)/n_y)./V_inf).^2;
            % Add all the pressure
            P = P_inf + Cp*q_inf; % total pressure

            % Determine color levels for stream function contours
            levmin = Psi(1,n_x); % defines the color levels -> trial and error to find a good representation
            levmax = Psi(n_y,n_x/50);
            levels = linspace(levmin,levmax,100)';

            %% Find the error associated with the 
            Psi_err = mean(abs((Psi_N-Psi)./Psi_N),'all');
            Phi_err = mean(abs((Phi_N-Phi)./Phi_N),'all');
            P_err = mean(abs((P_N-P)./P_N),'all');

            err = max([Psi_err Phi_err P_err]);
            err_vec(j) = err; 
            j = j + 1;

            if err > 0.005
                N = N + 100;
            end

        end
        

        N_final = N;

    % Plot
    figure
    contourf(x,y,Psi,75) % stream lines
    hold on 
    plot([0 c],[0 0],'k','linewidth',2) % airfoil
    axis equal
    axis(bounds)
    ylabel('y')
    xlabel('x')
    title(['Stream Lines for N = ' num2str(N_final) ' Vorticies at \alpha = 6^{o}']);

    % Plot equipotential lines at levels
    figure
    contourf(x,y,Phi,100) % equipotential lines
    hold on 
    plot([0 1.5],[0 0],'k','linewidth',2) % airfoil
    axis equal
    axis(bounds)
    ylabel('y')
    xlabel('x')
    title(['Equipotential Lines for N = ' num2str(N_final) ' Vorticies at \alpha = 6^{o}']);

    % Plot pressure contours
    figure
    contourf(x,y,P,100) % pressure contour
    hold on
    plot([0 c],[0 0],'k','linewidth',2) % airfoil
    axis equal
    axis(bounds)
    k = colorbar;
    k.Label.String = 'Pressure [Pa]';
    ylabel('y')
    xlabel('x')
    title(['Pressure Contours for N = ' num2str(N_final) ' Vorticies at \alpha = 6^{o}']);
        
        
%% If step == 2
        
    elseif step == 2
%       2 - regular function that is used to conducts study to determine 
%           changes in angle of attack
        
        % Domain Declaration
        x_min = -c/2;
        x_max = c+c/2;
        y_min = -c/2;
        y_max = c/2;
        bounds = [x_min x_max y_min y_max];
        N_final = N;

        % Define Number of Grid Points
        n_x = 100; % steps in the x direction
        n_y = 100; % steps in the y direction

        % Define Mesh
        [x,y] = meshgrid(linspace(x_min,x_max,n_x),linspace(y_min,y_max,n_y));

        % Define radius
        r =@(x1) sqrt((x-x1).^2 + (y).^2);

        % Define vortices
        delta_x = c/N;
        x_vortex = linspace(delta_x/2,c-delta_x,N);

        % Calculate vortex strength and circulation
        gamma = 2*alpha*V_inf*sqrt((1-(x_vortex/c))./(x_vortex/c)); % vortex strength
        Gamma = gamma*delta_x;                                      % circulation

        % Calculate psi for uniform flow (Eq. 3.55)
        Psi_uniform = V_inf*(y*cos(alpha)-x*sin(alpha));

        % Calculate phi for uniform flow (Eq. 3.53)
        Phi_uniform = V_inf*(x*cos(alpha)-y*sin(alpha));

        % Calculate psi for vortex (Eq. 3.114)
        Psi_vortex = 0;
        for i = 1:N
            % calculate the stream function of each vortices and add them
            % together to get the total stream function of the system
            Psi_vortex = Psi_vortex + (Gamma(i)*log(r(x_vortex(i))))/(2*pi);
        end

        % Calulate phi for vortex (Eq. 3.112)
        Phi_vortex = 0;
        for i = 1:N
            % calculate the velocity potential of each vortices and add them
            % together to get the total velocity potential of the system
            theta = atan2(-y,-x+x_vortex(i));
            Phi_vortex = Phi_vortex + -(Gamma(i)*theta)/(2*pi);
        end

        % Add all the stream functions
        Psi = Psi_vortex + Psi_uniform;

        % Add all the velocity potentials
        Phi = Phi_vortex + Phi_uniform;

        % Calculate P
        q_inf = (1/2)*rho_inf*V_inf^2; % dynamic pressure
        Cp = 1-(gradient(Phi,(x_max-x_min)/n_x,(y_max-y_min)/n_y)./V_inf).^2;
        % Add all the coeff. of pressure
        P = P_inf + Cp*q_inf; % total pressure
        
        % Determine color levels for stream function contours
        levmin = Psi(1,n_x); % defines the color levels -> trial and error to find a good representation
        levmax = Psi(n_y,n_x/2);
        levels = linspace(levmin,levmax,100)';
        err_vec = 1;
                %% If step == 3
    elseif step == 3
%       3 - Calculates and plots with inital N value input

        % Domain Declaration
        x_min = -c/2;
        x_max = c+c/2;
        y_min = -c/2;
        y_max = c/2;
        bounds = [x_min x_max y_min y_max];

        % Define Number of Grid Points
        n_x = 100; % steps in the x direction
        n_y = 100; % steps in the y direction

        % Define Mesh
        %Assume that 50000 is enough to be about 100% accurate
        N_2 = N;
        [x,y] = meshgrid(linspace(x_min,x_max,n_x),linspace(y_min,y_max,n_y));

        % Define radius
        r =@(x1) sqrt((x-x1).^2 + (y).^2);
        
        % Define Vorticies
        delta_x = c/N_2; %panel width
        x_vortex = linspace(delta_x/2,c-delta_x,N_2); %x-values along chord

        % Calculate vortex strength and circulation
        gamma = 2*alpha*V_inf*sqrt((1-(x_vortex/c))./(x_vortex/c)); % vortex strength
        Gamma = gamma*delta_x;                                      % circulation

        % Calculate psi for uniform flow (Eq. 3.55)
        Psi_uniform = V_inf*(y*cos(alpha)-x*sin(alpha));    

        % Calculate phi for uniform flow (Eq. 3.53)
        Phi_uniform = V_inf*(x*cos(alpha)-y*sin(alpha));    

        % Calculate psi for vortex (Eq. 3.114)
        Psi_vortex = 0;
        for i = 1:N_2
            % calculate the stream function of each vortex and add them
            % using super position to get the total stream function
            Psi_vortex = Psi_vortex + (Gamma(i)*log(r(x_vortex(i))))/(2*pi);
        end

        % Calulate phi for vortex (Eq. 3.112)
        Phi_vortex = 0;
        for i = 1:N_2
            % calculate the velocity potential of each vortices and add them
            % using super position to get the total velocity potential
            theta = atan2(-y,-x+x_vortex(i));
            Phi_vortex = Phi_vortex + -(Gamma(i)*theta)/(2*pi);
        end

        % Add all the stream functions
        Psi = Psi_vortex + Psi_uniform;

        % Add all the velocity potentials
        Phi = Phi_vortex + Phi_uniform; 

        % Calculate P
        q_inf = (1/2)*rho_inf*V_inf^2; % dynamic pressure
        Cp = 1-(gradient(Phi,(x_max-x_min)/n_x,(y_max-y_min)/n_y)./V_inf).^2; %  Cp = 1 - ([Vx Vy]/V)^2 from cylindrical flow
        % Add all the pressure
        P = P_inf + Cp*q_inf;
        
        % Determine color levels for stream function contours
        levmin = Psi(1,n_x); % defines the color levels -> trial and error to find a good representation
        levmax = Psi(n_y,n_x/50);
        levels = linspace(levmin,levmax,100)';
        
        N_final = N_2;
        
        err_vec = zeros(length(N:100:24000),1);
      
    % Plot
    figure
    contourf(x,y,Psi,75) % stream lines
    hold on 
    plot([0 c],[0 0],'k','linewidth',2) % airfoil
    axis equal
    axis(bounds)
    ylabel('y')
    xlabel('x')
    title(['Stream Lines for N = ' num2str(N_2) ' Vorticies at \alpha = 6^{o}']);

    % Plot equipotential lines at levels
    figure
    contourf(x,y,Phi,100) % equipotential lines
    hold on 
    plot([0 1.5],[0 0],'k','linewidth',2) % airfoil
    axis equal
    axis(bounds)
    ylabel('y')
    xlabel('x')
    title(['Equipotential Lines for N = ' num2str(N_2) ' Vorticies at \alpha = 6^{o}']);

    % Plot pressure contours
    figure
    contourf(x,y, P,100) % pressure contour
    hold on
    plot([0 c],[0 0],'k','linewidth',2) % airfoil
    axis equal
    axis(bounds)
    k = colorbar;
    k.Label.String = 'Pressure [Pa]';
    ylabel('y')
    xlabel('x')
    title(['Pressure Contours for N = ' num2str(N_2) ' Vorticies at \alpha = 6^{o}']);
        
    
        
end

