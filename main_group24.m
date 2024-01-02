clear
clc

addpath('.\data')
%addpath('.\scripts\')
%addpath('functions\')

stepsize = 0.1; % set lower value for a more refined smulation

% DATA

%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%

mu = 3.98600433e+5;  % Earth gravitational constant  [km^3/s^2]
c = 2.99792*10^8;  % Light velocity  [m/s] 

%%%%%%%%%%%%%%%%%%%%%%%%%% Orbit Initial Data %%%%%%%%%%%%%%%%%%%%%%%%%%

R_e = 6378;  % Earth's radius  [km]
a = 6939.063;  % Semi-major axis  [km]             
e = 0;  % Eccentricity  [-]             
i = deg2rad(97.5307);  % Inclination  [rad]           
theta0 = deg2rad(30);  % true anomaly initial condition  [rad] 
n = sqrt(mu/a^3);  % mean angular velocity  [rad/s]         
T = 2*pi/n;  % orbit period  [s]   

%%%%%%%%%%%%%%%%%%%%%% Spacecraft Characteristics %%%%%%%%%%%%%%%%%%%%%%

MB = [4; 0.1; 0.1; 0.3];% Main body  [kg; m; m; m] Mass, a, b, h
SP = [0.25; 0.1; 0.3; 0.01];% Solar panels  [kg; m; m; m] Mass, a, b, h     % (( ?? )) %  <--------  h = spessore?? serve?

x = 0.1; y = 0.1; z = 0.3; % Main body size [m]
mass = 4; % Main body mass  [kg]
% Inertia matrix main body 
Ix_mb = mass / 12*((y^2+z^2)); % [kg*m^2]
Iy_mb = mass / 12*((x^2+z^2));
Iz_mb = mass / 12*((y^2+x^2));
% Inertia matrix main body in s/c frame (Huygens-Steiner theorem)
Ix_hs_mb = Ix_mb + mass * ((z/2)/5)^2;
Iy_hs_mb = Iy_mb + mass * ((z/2)/5)^2;
Iz_hs_mb = Iz_mb;
% Define solar panel type A and B
x_spA = 0.3; x_spB = 0.1; % [m]
y_spA = 0.1; y_spB = 0.3; % [m]
z_sp = 0.002; % [m]
mass_sp = 0.25; % [kg]
% Inertia matrix panel A in panel frame
Ix_spA = mass_sp / 12*((y_spA^2+z_sp^2));
Iy_spA = mass_sp / 12*((x_spA^2+z_sp^2));
Iz_spA = mass_sp / 12*((y_spA^2+x_spA^2));
% Distance along z of the panels from center of mass
z_cg = 4/5*z/2;
% Inertia matrix panel A in s/c frame (Huygens-Steiner theorem)
Ix_hsA = Ix_spA + mass_sp * (z_cg)^2;
Iy_hsA = Iy_spA + mass_sp * ((z_cg)^2+(x/2+z/2)^2);
Iz_hsA = Iz_spA + mass_sp * (x/2+z/2)^2;
% Inertia matrix panel B in panel frame
Ix_spB = mass_sp / 12*((y_spB^2+z_sp^2));
Iy_spB = mass_sp / 12*((x_spB^2+z_sp^2));
Iz_spB = mass_sp / 12*((y_spB^2+x_spB^2));
% Inertia matrix panel B in s/c frame (Huygens-Steiner theorem)
Ix_hsB = Ix_spB + mass_sp * ((z_cg)^2+(x/2+z/2)^2);
Iy_hsB = Iy_spB + mass_sp * (z_cg)^2;
Iz_hsB = Iz_spB + mass_sp * (x/2+z/2)^2;
% Total inertia in body frame 
Ix = Ix_hs_mb + 2*Ix_hsA + 2*Ix_hsB; 
Iy = Iy_hs_mb + 2*Iy_hsA + 2*Iy_hsB;
Iz = Iz_hs_mb + 2*Iz_hsA + 2*Iz_hsB;
I_tot = [Ix 0 0; 0 Iy 0; 0 0 Iz];
invI_tot = inv(I_tot);     

%%%%%%%%%%%%%%%%%%%%%%% Sensors characteristics %%%%%%%%%%%%%%%%%%%%%%%

% Gyroscope STIM380H
GyroSampleRate = 262;  % Gyroscope sample rate  [Hz] 
ARW = 0.10;  % Angular Randon Walk gyroscope  [deg/sqrt(h)]
RRW = 0.5;  % Rate Randon Walk gyroscope  [deg/h]

% Magnetometer DTFM100S
MMAccuracy = 0.003;  % +-0.3%
MMSampleRate = 1;  % Magnetometer sample rate  [Hz] to speed up the simulation
MMnoise = 15e-9*[1; 1; 1];  % noise vector (bias)  [T]

% Extended State Observer
Lw = 0.8;
Ld = 1e-5;


%%%%%%%%%%%%%%%%%%%%%% Initial Conditions %%%%%%%%%%%%%%%%%%%%%%

wx0 = deg2rad(3);  %omega x  [rad/s]        
wy0 = deg2rad(11);  %omega y  [rad/s]                          
wz0 = deg2rad(14);  %omega z  [rad/s]                        
w0 = [wx0; wy0; wz0];  %omega vector  [rad/s] [3x1] 
A_BN0 = eye(3);
A0 = eye(3);
q0 = [0 0 0 1]'; % quaternions initial condition

%%%%%%%%%%%%%%%%%%%%%%%%%% Environment %%%%%%%%%%%%%%%%%%%%%%%%%%

% Magnetic Disturbance
% Roba necessaria per il Subsystem "Mag field senza Matlab fcn"
% n_E = 729.2e-07; % [rad/s]
% g1_0 = -29404.8 *1e-9;
% g1_1 = -1450.9 *1e-9;
% h1_1 = 4652.5 *1e-9;
% H_0 = ((g1_0^2)+(g1_1^2)+(h1_1^2))^0.5;

m = [0.1; 0.1; 0.1]; % [A*m^2] Residual  magnetic  induction  due  to  currents  in  the  satellite
[gn, gm, gvali, gsvi] = textread('igrfSg.txt','%f %f %f %f');
[hn, hm, hvali, hsvi] = textread('igrfSh.txt','%f %f %f %f');
G = [gn, gm, gvali, gsvi];
H = [hn, hm, hvali, hsvi];

% Solar Radiation Pressure
rho_d_MB = 0.1;
rho_s_MB = 0.5;
rho_s_SP = 0.1;
rho_d_SP = 0.1;


stopTime = T;




%% Variable Thrusters
bz = 0.1;                                                                   % braccio forze lungo asse z [m]
by = 0.1;                                                                   % braccio forze lungo asse y [m]

thrust.R = [bz -bz 0   0   0   0                                            % R matrix [3x6]
             0  0  bz -bz  0   0
             0  0  0   0  by -by];
thrust.R_pinv = pinv(thrust.R);                                             % pseudo invers R matrix
thrust.T_min = 10e-6;                                                     % minimum Thrust [N]
thrust.T_max = 500e-6;                                                    % maximum Thrust [N]
thrust.w = sum(null(thrust.R,'r'),2);                                       % non lo so 

%% Pointing control

k1 = 0.01;
k2 = 0.01;
