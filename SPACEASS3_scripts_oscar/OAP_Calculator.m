% AERO2705 Assignment 3 OAP Fuel Requirements | Author: Oscar Ansted
clear;
clc;
clf;

%% Global Parameters
deg = pi/180; % converts degrees to radians
mu = 3.986004418e+14; % standard gravitational parameter for Earth (m^3/s^2)
Re = 6378e+3; % radius of Earth (m)
longSyd = 152; % East longitude of Sydney (deg)
g0 = 9.81; % acceleration due to gravity @ SL (m/s/s)

%% Satellite to be serviced, spacecraft A
mA_E = 1300; % empty mass of spacecraft A (kg)

%% OAP, spacecraft C
% Input mass parameters
mC_P = 500; % propellant mass of OAP (kg)
mC_E = 50; % empty mass of OAP (kg) [including TCR and structural]
mC_0 = mC_P + mC_E; % net initial mass of the OAP (kg)

% Input specific impulse of chosen fuel 
Isp = 300; % specific impulse of OAP propellant (s)

% Input the required delta-v per year
delta_v_yr = 50; % delta-v needed to maintain the target per year (m/s)

% Input the number of years to support the target
nYears = 5; % number of years of support (years)

%% Calculate the mated spacecraft's parameters, at epoch (t = 0)
m0 = mA_E + mC_0; % initial empty mass of serviced spacecraft (kg)

%% Calculate the fuel required for the burns over each year


