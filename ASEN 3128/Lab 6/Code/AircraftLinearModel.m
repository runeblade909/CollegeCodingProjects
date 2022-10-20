function [Alon, Blon, Alat, Blat] = AircraftLinearModel(trim_definition, trim_variables, aircraft_parameters)
%
% STUDENT COMPLETE


u0 = trim_definition(1);
h0 = trim_definition(2);

alpha0 = trim_variables(1);
de0 = trim_variables(2);
dt0 = trim_variables(3);

theta0 = alpha0;

ap = aircraft_parameters;
rho = stdatmo(h0);



%%%%%%%%%%%%%%%%%%%%%%
%%% Longitudinal
%%%%%%%%%%%%%%%%%%%%%%

%%%% Trim values
CW0 = ap.W/((1/2)*rho*u0^2*ap.S);
CL0 = CW0*cos(theta0);
CD0 = ap.CDmin + ap.K*(CL0-ap.CLmin)^2;
CT0 = CD0 + CW0*sin(theta0);


%%%% Nondimensional stabiulity derivatives in body coordinates

%%%%% This is provided since we never discussed propulsion - Prof. Frew
dTdu = dt0*ap.Cprop*ap.Sprop*(ap.kmotor-2*u0+dt0*(-2*ap.kmotor+2*u0));
CXu = dTdu/(.5*rho*u0*ap.S)-2*CT0;

CDu = 0;%CDM*Ma;    % Compressibility only.  Ignore aeroelasticity and other effects (dynamic pressure and thrust).
CLu = 0;%CLM*Ma;
Cmu = 0;%CmM*Ma;

CZu = 0; %ASSUME

CZalpha = -CD0 - ap.CLalpha;
CXalpha = CL0*(1-2*ap.K*ap.CLalpha);

CZalphadot = -ap.CLalphadot;
CXalphadot = 0; %ASSUME

CZq = -ap.CLq;

CXq = 0; %ASSUME

% Longitudinal dimensional stability derivatives (from Etkin and Reid)
Xu = rho*u0*ap.S*CW0*sin(theta0) + 0.5*rho*u0*ap.S*CXu;
Zu = -rho*u0*ap.S*CW0*cos(theta0) + 0.5*rho*u0*ap.S*CZu;
Mu = 0.5*rho*u0*ap.S*ap.c*Cmu;


Xw = 0.5*rho*u0*ap.S*CXalpha;
Zw = 0.5*rho*u0*ap.S*CZalpha;
Mw = 0.5*rho*u0*ap.S*ap.c*ap.Cmalpha;

Xq = 0.25*rho*u0*ap.c*ap.S*CXq;
Zq = 0.25*rho*u0*ap.c*ap.S*CZq;
Mq = 0.25*rho*u0*ap.c^2*ap.S*ap.Cmq;

Xwdot = 0.25*rho*ap.c*ap.S*CXalphadot;
Zwdot = 0.25*rho*ap.c*ap.S*CZalphadot;
Mwdot = 0.25*rho*ap.c^2*ap.S*ap.Cmalphadot;

%Define Alon components outside matrix
i11 = Xu/ap.m;
i12 = Xw./ap.m;
i13 = 0;
i14 = -ap.g*cos(theta0);
i21 = Zu/(ap.m - Zwdot);
i22 = Zw/(ap.m-Zwdot);
i23 = (Zq + ap.m*u0)/(ap.m - Zwdot);
i24 = (-ap.m*ap.g*sin(theta0))/(ap.m - Zwdot);
i31 = 1/ap.Iy*(Mu + (Mwdot*Zu)/(ap.m - Zwdot));
i32 = 1/ap.Iy*(Mw + (Mwdot*Zw)/(ap.m - Zwdot));
i33 = 1/ap.Iy*(Mq + Mwdot*(Zq + ap.m*u0)/(ap.m - Zwdot));
i34 = (-Mwdot*ap.m*ap.g*sin(theta0))/(ap.Iy*(ap.m - Zwdot));
i41 = 0;
i42 = 0;
i43 = 1;
i44 = 0;

% Matrices
Alon = [i11 i12 i13 i14; i21 i22 i23 i24; i31 i32 i33 i34; i41 i42 i43 i44];

Blon = zeros(6,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Lateral
%%%%%%%%%%%%%%%%%%%%%%%%%%

% Lateral-directional dimensional stability derivatives
Yv = 0.5*rho*u0*ap.S*ap.CYbeta;
Yp = 0.25*rho*u0*ap.b*ap.S*ap.CYp;
Yr = 0.25*rho*u0*ap.b*ap.S*ap.CYr;

Lv = 0.5*rho*u0*ap.b*ap.S*ap.Clbeta;
Lp = 0.25*rho*u0*ap.b^2*ap.S*ap.Clp;
Lr = 0.25*rho*u0*ap.b^2*ap.S*ap.Clr;

Nv = 0.5*rho*u0*ap.b*ap.S*ap.Cnbeta;
Np = 0.25*rho*u0*ap.b^2*ap.S*ap.Cnp;
Nr = 0.25*rho*u0*ap.b^2*ap.S*ap.Cnr;

G = ap.Ix*ap.Iz-ap.Ixz^2;

G3=ap.Iz/G;
G4=ap.Ixz/G;
G8=ap.Ix/G;

j11 = Yv/ap.m;
j12 = Yp/ap.m;
j13 = Yr/ap.m - u0;
j14 = ap.g*cos(theta0);
j21 = G3*Lv + G4*Nv;
j22 = G3*Lp + G4*Np;
j23 = G3*Lr + G4*Nr;
j24 = 0;
j31 = G4*Lv + G8*Nv;
j32 = G4*Lp + G8*Np;
j33 = G4*Lr + G8*Nr;
j34 = 0;
j41 = 0;
j42 = 1;
j43 = tan(theta0);
j44 = 0;

Alat = [j11 j12 j13 j14; j21 j22 j23 j24; j31 j32 j33 j34; j41 j42 j43 j44];

Blat = zeros(6,2);


