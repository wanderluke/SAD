clear
clc

addpath('.\data')
%addpath('.\scripts\')
%addpath('functions\')

stepsize = 0.01; % set lower value for a more refined smulation

% DATA

%%%%%%%%%%%%%%%%%%%%%%%%%% Constants %%%%%%%%%%%%%%%%%%%%%%%%%%

mu = 3.98600433e+5;  % Earth gravitational constant  [km^3/s^2]
c = 2.99792*10^8;  % Light velocity  [m/s] 

%%%%%%%%%%%%%%%%%%%%%%%%%% Orbit Initial Data %%%%%%%%%%%%%%%%%%%%%%%%%%

R_e = 6378;  % Earth's radius  [km]
a = 6939.063;  % Semi-major axis  [km]             
e = 0.001304;  % Eccentricity  [-]             
i = deg2rad(97.5307);  % Inclination  [rad]           
theta0 = deg2rad(30);  % true anomaly initial condition  [rad] 
n = sqrt(mu/a^3);  % mean angular velocity  [rad/s]         
T = 2*pi/n;  % orbit period  [s]   

%%%%%%%%%%%%%%%%%%%%%% Spacecraft Characteristics %%%%%%%%%%%%%%%%%%%%%%

x = 0.1; y = 0.1; z = 0.3; % [m]
mass = 4; % [kg]
% Inertia body 
Ix_mb = mass / 12*((y^2+z^2)); % [kg*m^2]
Iy_mb = mass / 12*((x^2+z^2));
Iz_mb = mass / 12*((y^2+x^2));
% Define solar panel type A and B
x_spA = 0.3; x_spB = 0.1; % [m]
y_spA = 0.1; y_spB = 0.3; % [m]
z_sp = 0.002; % [m]
mass_sp = 0.25; % [kg]
% Inertia matrix panel A in panel frame
Ix_spA = mass_sp / 12*((y_spA^2+z_sp^2));
Iy_spA = mass_sp / 12*((x_spA^2+z_sp^2));
Iz_spA = mass_sp / 12*((y_spA^2+x_spA^2));
% Inertia matrix panel A in body frame (Huygens-Steiner theorem)
Ix_hsA = Ix_spA + mass_sp * (z/2)^2;
Iy_hsA = Iy_spA + mass_sp * ((z/2)^2+(x/2+z/2)^2);
Iz_hsA = Iz_spA + mass_sp * ((x/2+z/2)^2)^2;
% Inertia matrix panel B in panel frame
Ix_spB = mass_sp / 12*((y_spB^2+z_sp^2));
Iy_spB = mass_sp / 12*((x_spB^2+z_sp^2));
Iz_spB = mass_sp / 12*((y_spB^2+x_spB^2));
% Inertia matrix panel B in body frame (Huygens-Steiner theorem)
Ix_hsB = Ix_spB + mass_sp * ((z/2)^2+(x/2+z/2)^2);
Iy_hsB = Iy_spB + mass_sp * (z/2)^2;
Iz_hsB = Iz_spB + mass_sp * ((x/2+z/2)^2)^2;
% Total inertia in body frame 
Ix = Ix_mb + 2*Ix_hsA + 2*Ix_hsB; 
Iy = Iy_mb + 2*Iy_hsA + 2*Ix_hsB;
Iz = Iz_mb + 2*Iz_hsA + 2*Ix_hsB;
I_tot = [Ix 0 0; 0 Iy 0; 0 0 Iz];
invI_tot = inv(I_tot);     

%%%%%%%%%%%%%%%%%%%%%% Initial Conditions %%%%%%%%%%%%%%%%%%%%%%

wx0 = deg2rad(3);  %omega x  [rad/s]        
wy0 = deg2rad(11);  %omega y  [rad/s]                          
wz0 = deg2rad(14);  %omega z  [rad/s]                        
w0 = [wx0; wy0; wz0];  %omega vector  [rad/s] [3x1] 
A_BN0 = eye(3);

%%%%%%%%%%%%%%%%%%%%%%%%%% Environment %%%%%%%%%%%%%%%%%%%%%%%%%%

% Gravity Gradient
inertia1 = Iz-Iy;  % term to compute gravity gradient  [Kg*m^2]
inertia2 = Ix-Iz;  % term to compute gravity gradient  [Kg*m^2]
inertia3 = Iy-Ix;  % term to compute gravity gradient  [Kg*m^2]

% Magnetic Dsiturbance
m = [0.1; 0.1; 0.1];
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
