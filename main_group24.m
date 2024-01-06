clear
clc

%                                        /\   
%                                        |   z         
%                                        |                         
%                                    ________    
%                                   /        /  
%                                  /        / 
%                                 /        /                
%                                /        /      
%                               /        /
%                              /        /
%                             /        / 
%     +--------------------- +--------+ ---------------------+
%    /                      /|       / |                    /        y
%   /                      / |      /  |                   /   --------> 
%  +--------------------- +--------+ --|------------------+  
%                        /|  |    /|   |
%                       / |  |   / |   |
%                      /  |  |  /  |   |
%                     /   |  | /   |   |
%                    /    |  |/    |   |
%                   /     |  /     |   |
%                  /________/+-----|---+
%                         | /      |  /
%                         |/       | /
%                         +--------+
%                             /
%                            /  x
%                           \/
%



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
a = 6978;  % Semi-major axis  [km]             
e = 0;  % Eccentricity  [-]             
i = deg2rad(98);  % Inclination  [rad]           
theta0 = deg2rad(0);  % true anomaly initial condition  [rad] 
n = sqrt(mu/a^3);  % mean angular velocity  [rad/s]         
T = 2*pi/n;  % orbit period  [s]   

%%%%%%%%%%%%%%%%%%%%%% Spacecraft Characteristics %%%%%%%%%%%%%%%%%%%%%%

MB = [4; 0.1; 0.1; 0.3];% Main body  [kg; m; m; m] Mass, a, b, h
SP = [0.25; 0.1; 0.3; 0.01];% Solar panels  [kg; m; m; m] Mass, a, b, h  

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

%%%%%%%%%%%%%%%%%%%%%%% Surface of the main body %%%%%%%%%%%%%%%%%%%%%%% 
s1.n = [1; 0; 0];
s1.A = z*y;
s1.r = [x/2; 0; 0];
s2.n = [0; 1; 0];
s2.A = x*z;
s2.r = [0; y/2; 0];
s3.n = [0; 0; 1];
s3.A = x*y;
s3.r = [0; 0; 4/5*z/2]; 
s4.n = [-1; 0; 0];
s4.A = z*y;
s4.r = [-x/2; 0; 0];
s5.n = [0; -1; 0];
s5.A = x*z;
s5.r = [0; -y/2; 0];
s6.n = [0; 0; -1];
s6.A = x*y;
s6.r = [0; 0; -6/5*z/2];

pan.n1 = [0; 0; 1];
pan.n2 = [0; 0; -1];
pan.A = x_spA*y_spA;
pan.r = [0; 0; 4/5*z/2];

%%%%%%%%%%%%%%%%%%%%%%% Sensors characteristics %%%%%%%%%%%%%%%%%%%%%%%

% Gyroscope STIM380H
GyroSampleRate = 10;  % Gyroscope sample rate  [Hz] 
ARW = deg2rad(0.10);  % Angular Randon Walk gyroscope  [rad/s]
sigma_ARW = ARW/sqrt(3600);
sigma_n = sigma_ARW/sqrt(1/GyroSampleRate);
RRW = deg2rad(0.3);  % Rate Randon Walk gyroscope  [rad/s]
sigma_RRW = RRW/sqrt(3600);
sigma_b = sigma_RRW/sqrt(1/GyroSampleRate);
b0 = deg2rad(0.3)/sqrt(3600)*ones(3,1);

% Magnetometer DTFM100S
AccMagn = 0.5; 
MMSampleRate = 1;  % Magnetometer sample rate  [Hz] to speed up the simulation
MMnoise = 15e-9*[1; 1; 1];  % noise vector (bias)  [T]

% Sun Sensor
AccSunSens = 0.125;

%%%%%%%%%%%%%%%%%%%%%% Initial Conditions %%%%%%%%%%%%%%%%%%%%%%

wx0 = deg2rad(3);  %omega x  [rad/s]        
wy0 = deg2rad(11);  %omega y  [rad/s]                          
wz0 = deg2rad(14);  %omega z  [rad/s]                        
w0 = [wx0; wy0; wz0];  %omega vector  [rad/s] [3x1] 
q0 = [0 0 0 1]'; % quaternions initial condition

%%%%%%%%%%%%%%%%%%%%%%%%%% Environment %%%%%%%%%%%%%%%%%%%%%%%%%%

% Magnetic Disturbance
% Roba necessaria per il Subsystem "Mag field senza Matlab fcn"
n_E = 729.2e-07; % [rad/s]
g1_0 = -29404.8 *1e-9;
g1_1 = -1450.9 *1e-9;
h1_1 = 4652.5 *1e-9;
H_0 = ((g1_0^2)+(g1_1^2)+(h1_1^2))^0.5;

m = [0.1; 0.1; 0.1]; % [A*m^2] Residual  magnetic  induction  due  to  currents  in  the  satellite
[gn, gm, gvali, gsvi] = textread('igrfSg.txt','%f %f %f %f');
[hn, hm, hvali, hsvi] = textread('igrfSh.txt','%f %f %f %f');
G = [gn, gm, gvali, gsvi];
H = [hn, hm, hvali, hsvi];

% Solar Radiation Pressure
Fe = 1358; % [W/m^2] Power per unit surface
c = 299792458; % [m/s] Speed of light
P = Fe/c; % Average pressure due to radiation
rho_s_body = 0.5;
rho_d_body = 0.1*2/3; % Aggiungi commento
rho_s_sp = 0.1;
rho_d_sp = 0.1*2/3;

% Air Drag
rho_air = 6.967e-13; % [kg/m^3]
Cd = 2.2;
om_E = 0.000072921; % [rad/s]

stopTime = T;


%% Variable Thrusters
bz = z-(mass*z/2)/(mass+mass_sp);                                           % Force's arm along z-axis [m]
by = y/2;                                                                   % Force's arm along y-axis [m]

thrust.R = [bz -bz 0   0   0   0                                            % R matrix [3x6]
             0  0  bz -bz  0   0
             0  0  0   0  by -by];
thrust.R_pinv = pinv(thrust.R);                                             % pseudo invers R matrix
thrust.T_min = 10e-6;                                                       % minimum Thrust [N]
thrust.T_max = 0.35e-2;                                                     % maximum Thrust [N]
thrust.w = sum(null(thrust.R,'rational'),2);                                % null space vector
thrust.Gain=0.01;                                                            % controller gain
thrust.Filter_band=1e-3;                                                    % pre-filter omega band  
%% Pointing control

k1 = 5e-2;
k2 = 2.5e-2;

%% PLOTs
clc , close all
outDetumbling = sim('detumbling_group24.slx');

set(0,'defaultTextInterpreter','latex','defaultAxesFontSize',15);
set(0,'defaultAxesTickLabelInterpreter','latex');
set(0, 'defaultLegendInterpreter','latex');

figure('Name','Disturbances'),
hold on, grid on, box on
plot(outDetumbling.M_GG, 'linewidth',1.5);
plot(outDetumbling.M_SRP, 'linewidth',1.5);
plot(outDetumbling.M_MAG, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel(' [N m]')
legend('Gravity Gradient','SRP','Magnetique Torque')
xlim([0, outDetumbling.tout(end)])



figure('Name','Total Disturbances Torque'),
hold on, grid on, box on
plot(outDetumbling.disturbances, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel(' [N m]')
legend('$Td_x$','$Td_y$','$Td_z$')
xlim([0, outDetumbling.tout(end)])


figure('Name','Angular Velocity $\omega$'),
hold on, grid on, box on
plot(outDetumbling.w, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel(' [rad/s]')
legend('$\omega_x$','$\omega_y$','$\omega_z$')
xlim([0, outDetumbling.tout(end)])


figure('Name','Angular Velocity $\omega$ (detail)'),
hold on, grid on, box on
plot(outDetumbling.w, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel(' [rad/s]')
yline(1e-3,'k--'), yline(-1e-3,'k--')
legend('$\omega_x$','$\omega_y$','$\omega_z$','','')
xlim([0, outDetumbling.tout(end)])
ylim([-1e-2 1e-2])


figure('Name','Ideal Control Torque'),
hold on, grid on, box on
plot(outDetumbling.MC_ideal, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel(' [N m]')
legend('$MC_x$','$MC_y$','$MC_z$')
xlim([0, outDetumbling.tout(end)])


figure('Name','Ideal Control Torque (detail)'),
hold on, grid on, box on
plot(outDetumbling.MC_ideal, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel(' [N m]')
legend('$MC_x$','$MC_y$','$MC_z$')
xlim([1000, outDetumbling.tout(end)])


figure('Name','Thrust Level'),
hold on, grid on, box on
plot(outDetumbling.Thrust_Level, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel(' $[\%]$')
legend('$F_1$','$F_2$','$F_3$','$F_4$','$F_5$','$F_6$')
xlim([0 outDetumbling.tout(end)])


figure('Name','Thrust Level (detail)'),
hold on, grid on, box on
plot(outDetumbling.Thrust_Level, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel(' $[\%]$')
legend('$F_1$','$F_2$','$F_3$','$F_4$','$F_5$','$F_6$')
xlim([0 300])


%% Pointing plots

Pointing = sim("pointing_group24.slx");

figure('Name','Earth Pointing - Pointing Error'),
hold on, grid on, box on
plot(Pointing.point_err, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel('Error [deg]')

figure('Name','Tracking - Angular Velocity'),
hold on, grid on, box on
plot(Pointing.omega, 'linewidth',1.5);
xlabel('$t [s]$'), ylabel('$\omega$ [rad/s]')
legend('$\omega_x$','$\omega_y$','$\omega_z$')

figure('Name','Tracking - Mc'),
hold on, grid on, box on
plot(Pointing.Mc, 'linewidth',1.5);
legend('$M_{c,x}$',...
    '$M_{c,y}$',...
    '$M_{c,z}$')
xlabel('t $[s]$'), ylabel('$[N m]$')

figure('Name','Tracking - Thrust level'),
hold on, grid on, box on
plot(Pointing.th_level, 'linewidth',1.5);
legend('$th_{level}$')
xlabel('t $[s]$'), ylabel('$[%]$')
