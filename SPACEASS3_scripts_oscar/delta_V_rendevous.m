


function delta_v_in_gcef = delta_V_rendevous(coeA, coeB, nHours, nFig)

    % Global Parameters
    deg = pi/180; % converts degrees to radians
    mu = 3.986004418e+14; % standard gravitational parameter for Earth (m^3/s^2)

    % Calculate some additional parameters for the target orbit
    h_A = coeA(1); % angular momentum (m^2/s)
    e_A = coeA(2); 
    T_A = 2*pi/mu^2*(h_A/sqrt(1 - e_A^2))^3; % Target's period (s)
    n = 2*pi/T_A; % target's mean motion (rad/s)

    % GCEF State Vectors 
    [rA_0, vA_0] = coe2stateVecs(coeA); % initial state of the target
    [rB_0, vB_0] = coe2stateVecs(coeB); % initial state of the chaser

    % Obtain the initial position and velocity of B in the LVLH reference frame
    [dr0, dv0_neg, ~] = gcef_to_LVLH(rA_0, vA_0, rB_0, vB_0); % current state of the chaser in LVLH

    % Calculate the trajectory
    % Time of flight
    hour = 60*60; % one hour expressed in seconds (s) 
    t0 = 0; % initial time (s)
    tf = nHours*hour; % final time (s)
    t = tf - t0; % time of flight (s)

    % As docking is involved, the final relative velocity of the chaser must be
    % 0
    dvf_pos = 0; % final velocity of chaser (m/s)

    % Using the CW matrices, calculate the final position and velocity vectors
    % for the transfer trajectory in LVLH, dr- and dv-
    dv0_pos = -inv(Phi_rv(n, t))*Phi_rr(n, t)*dr0; % the velocity v(0+) after the first impulsive burn (m/s)

    % Now, calculate the final velocities
    dvf_neg = Phi_vr(n, t)*dr0 + Phi_vv(n, t)*dv0_pos; % final velocity of the trajectory before second impulsive burn (m/s)

    % Delta-v's required
    % Calculate the first burn's delta-v
    delta_v0 = dv0_pos - dv0_neg; % delta-v of burn 1 (m/s)

    % Calculate the second burn's delta-v
    delta_vf = dvf_pos - dvf_neg; % delta-v of burn 2 (m/s)

    % Calculate the total delta-v required
    delta_v_total = abs(delta_v0) + abs(delta_vf); % net delta-v needed for rendevous

    % Producing a simulation of the trajectory and rendevous
    % Calculate the CW matrices as functions of time for a set timespan
    dt = 10; % evaluate every dt seconds
    tspan = t0:dt:tf; % vector of times for the CW to be evaluated

    dx = zeros(1, length(tspan))';
    dy = zeros(1, length(tspan))';
    dz = zeros(1, length(tspan))';
    du = zeros(1, length(tspan))';
    dv = zeros(1, length(tspan))';
    dw = zeros(1, length(tspan))';
    drdot = zeros(1, length(tspan))';

    for count = 1:length(tspan)


       dr_t = Phi_rr(n, tspan(count))*(dr0) + Phi_rv(n, tspan(count))*dv0_pos; 

       dv_t = Phi_vr(n, tspan(count))*dr0 + Phi_vv(n, tspan(count))*dv0_pos;

       dx(count) = dr_t(1);
       dy(count) = dr_t(2);
       dz(count) = dr_t(3);

       du(count) = dv_t(1);
       dv(count) = dv_t(2);
       dw(count) = dv_t(3);

       drdot(count) = norm([dv_t(1) dv_t(2) dv_t(3)]);

    end



    % Plot the trajectory
    figure(nFig)

    plot3(dx, dy, dz, 'b');
    hold on

    % Plot the axes 
    XX = 5000;
    YY = 5000;
    ZZ = 5000;
    quiver3(0, 0, 0, XX, 0, 0, 'Color', 'k');
    quiver3(0, 0, 0, 0, YY, 0, 'Color', 'k');
    quiver3(0, 0, 0, 0, 0, ZZ, 'Color', 'k');
    text(XX, 0, 0, 'X');
    text(0, YY, 0, 'Y');
    text(0, 0, ZZ, 'Z');

    % Label the axes
    xlabel('x_{LVLH} (m)');
    ylabel('y_{LVLH} (m)');
    zlabel('z_{LVLH} (m)');

    % Label the initial position of the chaser
    text(dx(1), dy(1), dz(1), 'x Start');
    axis on 
    grid on
    grid minor
    axis equal

    % Convert delta-v's to GCEF
    [~, DV_gcef] = LVLH_to_gcef(coeA, dr0, dv0_pos);
    [~, DV_gcef2] = LVLH_to_gcef(coeA, dr0, dv0_neg);

    delta_v_in_gcef = DV_gcef - DV_gcef2;




end