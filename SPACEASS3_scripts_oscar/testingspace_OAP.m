clear;
clc;
clf;
% Global Parameters
deg = pi/180; % converts degrees to radians
mu = 3.986004418e+14; % standard gravitational parameter for Earth (m^3/s^2)
Re = 6378e+3; % radius of Earth (m)
longSyd = 152; % East longitude of Sydney (deg)


% Initialise orbital parameters of Target, A
h_A = 1.296404445843445e+11; % (m^2/s)
e_A = 0.0;
i_A = 0*deg;
Omega_A = 0*deg;
omega_A = longSyd*deg;
theta_A = 0*deg;

coeA = [h_A, e_A, i_A, Omega_A, omega_A, theta_A]; 

% Calculate some additional parameters for the target orbit
T_A = 2*pi/mu^2*(h_A/sqrt(1 - e_A^2))^3; % Target's period (s)
n = 2*pi/T_A; % target's mean motion (rad/s)

% Initialise orbital parameters of Chaser, B 
h_B = 1.296404445843445e+11;
e_B = 0.00045;
i_B = 0*deg;
Omega_B = 0*deg;
omega_B = longSyd*deg;
theta_B = 0*deg;

coeB = [h_B, e_B, i_B, Omega_B, omega_B, theta_B]; 

% Convert COE's to RV's
[rA, vA] = coe2stateVecs(coeA);
[rB, vB] = coe2stateVecs(coeB);
disp(rB);
disp(vB);

% Convert to LVLH
[rrelB, vrelB, ~] = gcef_to_LVLH(rA, vA, rB, vB);

% Convert back to GCEF
[rB2, vB2] = LVLH_to_gcef(coeA, rrelB, vrelB);

disp(rB2);
disp(vB2);


