% Calculate the injection orbit
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

%% Initial RSV orbit (Orbit 1)
i1 = 27*deg; % inclination of orbit 1 (rad)
r_p1 = 250*10^3 + R_e; % apogee radius (m)
r_a1 = 14880*10^3 + R_e; % perigee radius (m)

% Adjust these elements
Omega1 = 0*deg; % ra of asc node (rad)
omega1 = 0*deg; % arg of perigee (rad)
theta1 = 0*deg; % true anomaly at the burn (rad)

% Calculate the additional elements
a1 = (r_a1 + r_p1)/2; % semi-major axis (m)
e1 = r_a1/a1 - 1; % eccentricity 
h1 =  sqrt(mu*a1*(1 - e1^2)); % specific angular momentum of geo orbit (m^2/s)
T1 = sqrt(a1^3*(4*pi^2)/mu); % period (s)
%% ODE this orbit over one period for visualisation
t0 = 0; 
tf = 1*T1;
tspan = [t0, tf]; % timespan vector (s)

% Convert to an initial state vector
coe = [h1, e1, Omega1, i1, omega1, theta1];
[R0, V0] = coe2stateVecs(coe);
f0 = [R0, V0]; % initial state of the RSV

% Call the ODE
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
[t, f] = ode45(@(t, f)twoBody(t, f), tspan, f0, options);

% Plot the orbit
figure(1)
earth_sphere('km'); 
hold on
plot3(f(:, 1)*km, f(:, 2)*km, f(:, 3)*km, 'r');
axis equal
xlabel('X (km)');
ylabel('Y (km)');
zlabel('Z (km)');
hold on
%% Target geostationary orbit (Orbit 2) [do we assume a phaser is needed?]
T2 = 86164; % period of one sidereal day (s)
a2 = (T2^2*mu/(4*pi^2))^(1/3); % semi-major axis (m)
e2 = 0; % eccentricity
h2 = sqrt(mu*a2*(1 - e2^2)); % specific angular momentum of geo orbit (m^2/s)
i2 = 0*deg; % inclination of GEO (rad)
Omega2 = 0*deg; % right ascension is undefined as i = 0
omega2 = longSyd*deg; % define the argument of perigee (rad) 
theta2 = 0*deg; % true anomaly (rad)

%% Plot the geostationary orbit on the same axes
t0 = 0; 
tf = 1*T2;
tspan = [t0, tf]; % timespan vector (s)

% Convert to an initial state vector
coe = [h2, e2, Omega2, i2, omega2, theta2];
[R0, V0] = coe2stateVecs(coe);
f0 = [R0, V0]; % initial state of the RSV

% Call the ODE
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
[t, f] = ode45(@(t, f)twoBody(t, f), tspan, f0, options);

% plot on figure 1
plot3(f(:, 1)*km, f(:, 2)*km, f(:, 3)*km, 'b');
hold on % wait to plot the hohmann transfer
%% Trajectory Calculations
% Hohmann is more efficient, so firstly transfer to geostationary altitude
% in the 27 deg plane

apogeeBurn = false;
[delta_vs, hohmann] = coplanarHohmann(a1, e1, a2, e2, 0, 0, apogeeBurn);

% Calculate the Hohmann transfer orbit and propagate in figure 1
h = hohmann(1);
e = hohmann(2);
a = hohmann(1)^2/(mu*(1 - hohmann(2)^2)); % semi-major axis (m)
T = sqrt(a^3*(4*pi^2)/mu); % period (s)

coe_transfer = [h, e, Omega1, i1, omega1, theta1];
[R_hohmann, V_hohmann] = coe2stateVecs(coe_transfer);


%% ODE the Hohmann transfer
t0 = 0; 
tf = 1*T/2;
tspan = [t0, tf]; % timespan vector (s)

% Convert to an initial state vector
f0 = [R_hohmann, V_hohmann]; % initial state of the RSV

% Call the ODE
options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
[t, f] = ode45(@(t, f)twoBody(t, f), tspan, f0, options);

% plot on figure 1
plot3(f(:, 1)*km, f(:, 2)*km, f(:, 3)*km, 'g');

%% Calculate the plane change delta-v
delta = 27*deg; % change in plane angle
v_geo = h2/a2; % target velocity
delta_v_planechange = commonApsePlaneChange(v_geo, v_geo, delta); % delta-v required for plane change

%% As time of epoch is now specified, assume a 180 degree phaser is required to meet sydney
% However, in the presentation, we will say the epoch is calibrated to
% ensure there is no phasing needed.

nPeriods = 6; % define the number of periods the phasing orbit completes to synchronise
dLong = 180*deg; % change in longitude (rad)
T_geo = 86164; % geostationary orbit period
omega_geo = 2*pi/T_geo; % angular speed of geo (rad/s)

% Calculate the elements of the phasing orbit
T_phase = (nPeriods*2*pi + dLong)/(omega_geo*nPeriods); % required period of phasing orbit
a_phase = ((mu*T_phase^2)/(4*pi^2))^(1/3); % semi-major axis (m)

r_p_phase = a2; % the perigee of the phaser is the radius of the circular geostationary orbit (m)
r_a_phase = 2*a_phase - r_p_phase; % apogee of transfer orbit (m)
e_phase = (r_a_phase - r_p_phase)/(r_a_phase + r_p_phase); % eccentricity of phaser 

h_phase = sqrt(mu*a_phase*(1 - e_phase^2)); % angular momentum of phaser 

% Now calculate the velocities required to compute the delta-v
v_B_2 = v_geo; % the initial orbit, orbit 2, has geostationary velocity (m/s)
v_B_phase = h_phase/r_p_phase; % velocity required for phaser (m/s)

% Calculate the delta-v's, and the net delta-v for phasing
delta_v_phase_1 = v_B_phase - v_B_2; % delta-v for first burn (m/s)
delta_v_phase_2 = -(v_B_phase - v_B_2); % delta-v for first burn (m/s)

delta_v_phase_total = abs(delta_v_phase_1) + abs(delta_v_phase_2); % net delta-v for phaser

% finally, calculate precession rate relative to the ground station
precessionLongitude = (dLong/(nPeriods*T_phase))/deg;

%% Print the results here
fprintf('INJECTION TRAJECTORY\n\n');
fprintf('BURN 1: From injection orbit to %2.0f deg inclined geosynchronous orbit\n', i1/deg);
fprintf('BURN 1.1\n');
fprintf('Delta-v: %f m/s\n', delta_vs(1));
fprintf('BURN 1.2\n');
fprintf('Delta-v: %f m/s\n', delta_vs(2));
fprintf('Net total delta-v: %f m/s\n', sum(delta_vs));
fprintf('Time elapsed: %f hours\n\n', T/60/60/2);

if omega1 == 0 || omega1 == pi
    
    fprintf('As position B lies at coplanar intersect, induce BURN 2 simultaneous to BURN 1\n');
    waitTime = 0;
else
    % Calculate the wait time
    waitTime = (pi - omega1)/omega_geo; % wait time (s)
    % Convert to hours
    waitTime = waitTime/60/60; % wait time (hrs)
    
    % display the wait time
    fprintf('Wait %4.2f hours after position B is reached to engage BURN 2\n\n', waitTime);
    
end

% Print onwards
fprintf('BURN 2: Plane change to geostationary orbit\n');
fprintf('Delta-v: %f m/s, at %f degrees to the plane of orbit 2\n\n', delta_v_planechange, i1/deg);
fprintf('RUNNING TOTALS // FINAL TOTALS IF NO PHASING IS REQUIRED\n');
fprintf('Delta_v: %f m/s\n', sum(delta_vs) + delta_v_planechange);
fprintf('Time elapsed: %f hours\n\n', waitTime + T/60/60/2);
fprintf('Proceed with BURN 3 if out of phase by specified time\n\n');
fprintf('BURN 3: Phasing orbit by %f degrees\n', dLong/deg);
fprintf('BURN 3.1\n');
fprintf('Delta-v: %f m/s\n', delta_v_phase_1);
fprintf('Delta-v: %f m/s\n', delta_v_phase_2);
fprintf('Net total delta-v: %f m/s\n', delta_v_phase_total);
fprintf('Time elapsed: %f hours\n\n', nPeriods*T_phase/60/60);
fprintf('FINAL TOTALS:\n');
fprintf('DELTA-V: %f m/s\n', sum(delta_vs) + delta_v_phase_total + delta_v_planechange);
fprintf('TIME ELAPSED: %f hours\n', nPeriods*T_phase/60/60 + waitTime +  T/60/60/2);

