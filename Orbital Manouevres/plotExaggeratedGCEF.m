

function plotExaggeratedGCEF(tspan, coeA, drB_0, dvB_0, R0, V0, exg, nFig, color)

    % Orbit A is the target geostationary
    % Orbit B is the chaser orbit
    
    drB_0_exg = drB_0*exg; % initial position vector in LVLH (m)
    dvB_0_exg = dvB_0*exg; % initial velocity vector in LVLH (m/s)
    
    [rB_0, vB_0] = LVLH_to_gcef(coeA, drB_0_exg', dvB_0_exg');
    
    f0_A = [R0, V0];
    f0_B = [rB_0, vB_0]; % compile the initial state vector for chaser

    options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
    [~, fA] = ode45(@(t, f)twoBody(t, f), tspan, f0_A, options);
    [~, fB] = ode45(@(t, f)twoBody(t, f), tspan, f0_B, options);

    % Plot the desired orbit in nFig
    figure(nFig)
    
    % Plot the earth
    earth_sphere('km');
    hold on
    km = 10^-3; % converts metres to kilometres
    plot3(fA(:, 1)*km, fA(:, 2)*km, fA(:, 3)*km, 'r');
    hold on
    plot3(fB(:, 1)*km, fB(:, 2)*km, fB(:, 3)*km, color);
    hold on
    text(fB(1, 1)*km, fB(1, 2)*km, fB(1, 3)*km, 'o Start');
    hold on
    quiver3(0, 0, 0, 40000, 0, 0, 'k');
    legend('', 'Target orbit', 'Chaser orbit');

    grid on
    grid minor
    axis equal

end