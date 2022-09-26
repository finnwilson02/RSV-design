% AERO2705 Assignment 3 Fly-By Calculator | Author: Oscar Ansted
clear;
clc;
clf; 

%% Global Parameters
deg = pi/180; % converts degrees to radians
mu = 3.986004418e+14; % standard gravitational parameter for Earth (m^3/s^2)
Re = 6378e+3; % radius of Earth (m)
longSyd = 152; % East longitude of Sydney (deg)


%% Initialise orbital parameters of Target, A
h_A = 1.296404445843445e+11; % (m^2/s)
e_A = 0.0;
i_A = 0*deg;
Omega_A = 0*deg;
omega_A = longSyd*deg;
theta_A = 0*deg;

coeA = [h_A, e_A, i_A, Omega_A, omega_A, theta_A]; 

%% Initialise orbital parameters of Chaser, B
h_B = 1.296404445843445e+11;
e_B = 0.00005;
i_B = 0*deg;
Omega_B = 0*deg;
omega_B = longSyd*deg;
theta_B = 0*deg;

coeB = [h_B, e_B, i_B, Omega_B, omega_B, theta_B]; 

%% GCEF State Vectors 
[rA_0, vA_0] = coe2stateVecs(coeA); % initial state of the target
[rB_0, vB_0] = coe2stateVecs(coeB); % initial state of the chaser

%% Initialise the orbit simulation in LVLH
% Set the initial time
t0 = 0; % initial time (s)
T_A = 2*pi/mu^2*(h_A/sqrt(1 - e_A^2))^3; % Target's period (s)

% Plot over m periods
m = 1; % Plot over m orbits 

% Specify the number of time steps per period of the target's orbit
n = 1000; 

% Calculate the timestep size
dt = T_A/n; % timestep (s)

% Initialise the time so that for the first increment of the loop, t = 0
t = -dt; % initial time (s)

% Initialise empty vectors to store the data needed for the plots
x = zeros(1, m*n);
y = zeros(1, m*n);
z = zeros(1, m*n);
r = zeros(1, m*n);
T = zeros(1, m*n);
% Iteratively evaluate the position and velocity vectors to produce a plot
% of the orbit
for count = 1:m*n
    
    % Calculate the time at this step
    t = t + dt; % current time (s)
    
    % Calculate the new position and velocity for both orbits in GCEF
    [rA, vA] = rv_from_r0v0(rA_0, vA_0, t); % state of the target
    [rB, vB] = rv_from_r0v0(rB_0, vB_0, t); % state of the chaser
    
    % Calculate the position, velocity and acceleration of the chaser
    % relative to the target in LVLH frame
    [rrel, vrel, arel] = gcef_to_LVLH(rA, vA, rB, vB); % current state of the chaser in LVLH
    
    % Store the current values in respective vectors to plot
    x(count) = rrel(1); % X-pos in LVLH (m)
    y(count) = rrel(2); % Y-pos in LVLH (m)
    z(count) = rrel(3); % Z-pos in LVLH (m)
    r(count) = norm(rrel); % distance between bodies (m)
    T(count) = t; % time (s)
     
end
    
    
%% Plot the results
figure(1)
plot3(x, y, z, 'b');

hold on
axis equal
axis on
grid on 
axis equal

% km = 1e-3; % metres to km 
% plot3(x, y, z, 'b');
% xlabel('y_{LVLH} (m)');
% ylabel('x_{LVLH} (m)');
% zlabel('z_{LVLH} (m)');
% 
% hold on
% % Plot the coordinate axes of the LVLH frame
% XX = 1.3*max(x); % (m)
% YY = 1.3*max(y); % (m)
% quiver3(0, 0, 0, XX, 0, 0, 'Color', 'k');
% quiver3(0, 0, 0, 0, YY, 0, 'Color', 'k');
% text(XX, 0, 0, 'X');
% text(0, YY, 0, 'Y');
% 
% if sum(z) == 0 
%     ZZ = 0; % (m)
% else
%     ZZ = 1.3*max(f(:, 3)); % (m)
%     quiver3(0, 0, 0, 0, 0, ZZ, 'Color', 'k');
%     text(0, 0, ZZ, 'Z');
% end


text(x(1), y(1), z(1), 'x Start');

% Plot the coordinate axes of the LVLH frame
XX = 1.3*max(x); % (m)
YY = 1.3*max(y); % (m)
quiver3(0, 0, 0, XX, 0, 0, 'Color', 'k');
quiver3(0, 0, 0, 0, YY, 0, 'Color', 'k');
text(XX, 0, 0, 'Y');
text(0, YY, 0, 'X');

if sum(z) == 0 
    ZZ = 0; % (km)
else
    ZZ = 1.3*max(z); % (km)
    quiver3(0, 0, 0, 0, 0, ZZ, 'Color', 'k');
    text(0, 0, ZZ, 'Z');
end
disp('closest pass (m)');
disp(min(r));

