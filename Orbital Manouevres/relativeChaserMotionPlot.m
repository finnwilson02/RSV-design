% This function plots the relative motion of a chaser (B) to a target (A)
% in LVLH reference frame, and also produces the relative velocity plot for
% each point

function [t, fFB1] = relativeChaserMotionPlot(tspan, drFB1_0, dvFB1_0, RA_0, VA_0, nFigs, color)

    f0FB1 = [drFB1_0 dvFB1_0]'; % initial state vector for numerical integration

    % Apply the ODE45 solver
    options = odeset('RelTol', 1e-10, 'AbsTol', 1e-10);
    [t, fFB1] = ode45(@(t, f)linearRelMotionRates(t, f, RA_0, VA_0), tspan, f0FB1, options); 

    % Obtain the relative velocity vector for motion blur calculations
    relVel_FB1 = fFB1(:, 4:6); % vector of all velocities at time throughout the motion

    % Loop to calculate the normalised velocities at each sample point
    vrel = zeros(length(relVel_FB1), 1); 
    for k = 1:length(relVel_FB1)
        vrel(k) = norm(relVel_FB1(k, :));
    end

    % Plot the output 
    figure(nFigs(1)) % Figure 1, the position of the chaser relative to the target
    km = 1e-3; % metres to km 
    plot3(fFB1(:, 2)*km, fFB1(:, 1)*km, fFB1(:, 3)*km, color);
    text(fFB1(1, 2)*km, fFB1(1, 1)*km, fFB1(1, 3)*km, 'x Start');
    xlabel('y_{LVLH} (km)');
    ylabel('x_{LVLH} (km)');
    zlabel('z_{LVLH} (km)');

    hold on
    % Plot the coordinate axes of the LVLH frame
    XX = 1.3*max(fFB1(:,1)*km); % (km)
    YY = 1.3*max(fFB1(:, 2)*km); % (km)
    quiver3(0, 0, 0, 0, XX, 0, 'Color', 'k');
    quiver3(0, 0, 0, YY, 0, 0, 'Color', 'k');
    text(YY, 0, 0, 'Y');
    text(0, XX, 0, 'X');

    if sum(fFB1(:, 3)) == 0 
        ZZ = 0; % (km)
    else
        ZZ = 1.3*max(fFB1(:, 3)*km); % (km)
        quiver3(0, 0, 0, 0, 0, ZZ, 'Color', 'k');
        text(0, 0, ZZ, 'Z');
    end
    axis equal
    grid on 
    grid minor

    figure(nFigs(2)) % figure 2, the speed of the chaser relative to the target over time
    plot(t, vrel, 'b');
    xlabel('time (s)');
    ylabel('v_{rel} (m/s)');

    grid on 
    grid minor


end