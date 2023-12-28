

clear,clc


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




