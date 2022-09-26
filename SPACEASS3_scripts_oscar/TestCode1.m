% AERO2705 Assignment 3 Linearlised Equations of relative motion in orbit
% Calculator | Author: Oscar Ansted
clear;
clc;
clf;
%% Global Parameters
mu = 3.986004418e+14; % standard gravitational parameter for Earth (m^3/s^2)
Re = 6378e+3; % radius of Earth (m)
deg = pi/180; % converts degrees to radians

%% Target vehicle input parameters
rp = Re + 300e+3; % perigee radius of target vehicle, A (m)
e = 0.0; % target eccentricity
i = 15*deg; % target inclination (rad)
Omega = 0*deg; % target ra of asc node (rad)
omega = 0*deg; % target arg of perigee (rad)
theta = 0*deg; % target true anomaly (rad)

% Calculate the additional orbital elements of the target
ra = rp*(1 + e)/(1 - e); % radius of apogee of target (m)
a = (ra + rp)/2; % semi-major axis of target (m)
T = sqrt((4*pi^2/mu)*a^3); % period of target orbit (s)
n = 2*pi/T; % mean motion of target orbit (rad/s)
h = sqrt(2*mu*rp*ra/(ra + rp)); % specific angular momentum of target (m^2/s)

% Compile the COE vector
coe = [h e Omega i omega theta]; 

%% Chaser vehicle input parameters & state vectors
dr_0 = [1000, 0, 0]; % initial position vector in LVLH (m)
dv_0 = [5, 0, 5]; % initial velocity vector in LVLH (m/s)

% Calculate the target's initial state vectors using the COE's
[R0, V0] = coe2stateVecs(coe); % extract the state vectors of the target 

%% Numerical integration using Runge-Kutta method
% Construct the timespan vector
t0 = 0; % initial time (s)
tf = 5*T; % final time (s)
tspan = [t0, tf]; % timespan (s)

f0 = [dr_0 dv_0 R0 V0]'; % initial state vector for numerical integration

% Apply the ODE45 solver
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
[t, f] = ode45(@(t, f)linearRelMotionRatesTEST(t, f), tspan, f0, options); 

% Obtain the relative velocity vector for motion blur calculations
VREL = f(:, 4:6); % vector of all velocities at time throughout the motion
vrel = zeros(length(VREL), 1); % vector of the norm velocity at each time
for k = 1:length(VREL)
    vrel(k) = norm(VREL(k, :));
end


%% Plot the output 
figure(1) % Figure 1, the position of the chaser relative to the target
km = 1e-3; % metres to km 
plot3(f(:, 2)*km, f(:, 1)*km, f(:, 3)*km, 'b');
text(f(1, 2)*km, f(1, 1)*km, f(1, 3)*km, 'x Start');
xlabel('y_{LVLH} (km)');
ylabel('x_{LVLH} (km)');
zlabel('z_{LVLH} (km)');

grid on 
grid minor

% figure(2) % figure 2, the velocity of the chaser relative to the target along the LVLH axes
% plot3(f(:, 5)*km, f(:, 4)*km, f(:, 6)*km, 'b');
% text(f(1, 5)*km, f(1, 4)*km, f(1, 6)*km, 'x Start');
% xlabel('v_y (km/s)');
% ylabel('v_x (km/s)');
% zlabel('v_{z} (km/s)');
% 
% grid on 
% grid minor 

figure(3) % figure 3, the speed of the chaser relative to the target over time
plot(t, vrel, 'b');
xlabel('time (s)');
ylabel('v_{rel} (m/s)');

grid on 
grid minor