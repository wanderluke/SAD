
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

MB = [4; 0.1; 0.1; 0.3];% Main body  [kg; m; m; m] Mass, a, b, h
SP = [1; 0.3; 0.1; 0.01];% Solar panels  [kg; m; m; m] Mass, a, b, h
% Inertia Moments main body
Ix_mb = MB(1)/12 * (MB(3)^2 + MB(2)^2); %  [Kg*m^2]
Iy_mb = MB(1)/12 * (MB(3)^2 + MB(4)^2); %  [Kg*m^2]
Iz_mb = MB(1)/12 * (MB(2)^2 + MB(4)^2); %  [Kg*m^2]

% Inertia Moments solar panels
Ix_sp = (2*SP(1))/12 * (SP(3)^2 + SP(2)^2); %  [Kg*m^2]
Iy_sp = (2*SP(1))/12 * (SP(3)^2 + SP(4)^2); %  [Kg*m^2]
Iz_sp = (2*SP(1))/12 * (SP(2)^2 + SP(4)^2); %  [Kg*m^2]

% Terms for huygens steiner theorem
Ix_hs = Ix_sp + 2*SP(1)*(2.5)^2; %  [Kg*m^2]
Iy_hs = Iy_sp; % distance = 0  [Kg*m^2]
Iz_hs = Iz_sp + 2*SP(1)*(2.5)^2; %  [Kg*m^2]

% total inertia moments
Ix = 57749.855*10^-6;
Iy = 41435.420*10^-6;
Iz = 24624.309*10^-6;

Inertia_matrix = [Ix 0 0; 0 Iy 0; 0 0 Iz];  % inertia matrix  [Kg*m^2]
Inverse_inertia = inv(Inertia_matrix);  % inverse inertia matrix  [Kg*m^2]        

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