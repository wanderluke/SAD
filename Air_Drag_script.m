%% Air_Drag
n1 = [1 0 0]';      n4 = -n1;                                               % Normals to the surfaces in Body RF
n2 = [0 1 0]';      n5 = -n2;
n3 = [0 0 1]';      n6 = -n3;
                                                                            
n7 = [1 0 0]';      n8 = -n7;                                               % Solar panels
n9 = [1 0 0]';      n10 =-n9;

%% shift mass center wrt geometric center

r1 =[10 0 0]'*1e-2;      r4 =-r1;                                           % shift between CM - CG [m^2]
r2 =[0,10,0]'*1e-2;      r5 =-r2;
r3 =[0,0,10]'*1e-2;      r6 =-r3;

r7 =[0,0,45]'*1e-2;      r8 =r7;
r9 =[0,0,-45]'*1e-2;     r10 =r9;

%% area

A1 = 6e-2; A4 = 6e-2;                                                       % Lateral faces area [m^2]
A2 = 6e-2; A5 = 6e-2; 

A3 = 4e-2; A6 = 4e-2;                                                       % Top faces area [m^2]

A7 = 12e-2;    A8 = 12e-2;                                                  % Solar panels area [m^2]                                  
A9 = 12e-2;    A10 = 12e-2; 

%% global data
J = diag([100.9 25.1 91.6])*1e-2;                                           % Inertia matrix in body frame [Kg m^2]
Ix = J(1,1); Iy = J(2,2); Iz = J(3,3);

rho_air = 1.65e-13;                                                         %  air density constant [Kg/m^3]

cd = 2.2;                                                                   % Drag coefficient [-]

w_earth =[0 0 7.2921e-5]';                                                  % Earth rotation around its axis [rad/s]

% sat orbital parameters
a = 6939;                                                                   % semi-major axis [Km]
e = 0;                                                                      % eccentricity [-]
i = deg2rad(97.5307);                                                       % orbit inclination [rad]
n = 0.0011;                                                                 % mean orbit angular velocity [rad/s]
muE = astroConstants(13);
T_1orbit =2*pi * sqrt(a^3/muE);                                            % Period of 1 orbit [s]

w0 = [0 0 0]';                                                              % initial sat angular velocity in body frame [rad/s]
q0 = [ 1 1 1 1 ]';
r_earth = astroConstants(2);                                                % 1 AU [Km]
epsilon = deg2rad(23.45);                                                   % ecliptic inclination [rad] 

A0 = eye(3);











