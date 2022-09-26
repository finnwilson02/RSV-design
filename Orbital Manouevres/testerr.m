
clear;
clc;
clf;
%% Global parameters
R_e = 6.3781e+6; % Mean radius of Earth (m)
mu = 3.986004418e+14; % standard gravitational parameter for Earth (m^3/s^2)
J2 = 0.00108263; % second zonal harmonics effect 
deg = pi/180; % converts degrees to radians
km = 1e-3; % converts metres to km
longSyd = 152; % longitude of sydney (deg)

%% orbit 1
r_p1 = 480*10^3 + R_e; 
r_a1 = 800*10^3 + R_e;
a1 = (r_p1 + r_a1)/2; 
e1 = abs((r_a1 - r_p1)/(r_a1+r_p1)); 

r_p2 = 16000*10^3 + R_e; 
r_a2 = 16000*10^3 + R_e;
a2 = (r_p2 + r_a2)/2; 
e2 = abs((r_a2 - r_p2)/(r_a2 + r_p2)); 


delta_vs = coplanarHohmann(a1, e1, a2, e2, 0, 0, 0);

