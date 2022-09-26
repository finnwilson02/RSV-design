% Function: rates
% Rates function for ODE45 numerical integration of the linearised
% equations of relative motion for orbiting bodies


function dydt = linearRelMotionRates(t, y, R0, V0)

    % Global parameters
    mu = 3.986004418e+14; % standard gravitational parameter for Earth (m^3/s^2)

    % Extract the input vector's variables
    dx = y(1); % x-pos of chaser (m)
    dy = y(2); % y-pos of chaser (m)
    dz = y(3); % z-pos of chaser (m)
    dvx = y(4); % x-vel of chaser (m)
    dvy = y(5); % x-vel of chaser (m)
    dvz = y(6); % x-vel of chaser (m)
    
    
    % Calculate the updated state vector using the universal anomaly
    % method
    [R, V] = rv_from_r0v0(R0, V0, t);
    
    % Extract the coordinates from the current time's state vectors
    X = R(1); % target X-coordinate in gcef (m)
    Y = R(2); % target Y-coordinate in gcef (m)
    Z = R(3); % target Z-coordinate in gcef (m)

    VX = V(1); % target X-velocity in gcef (m)
    VY = V(2); % target Y-velocity in gcef (m)
    VZ = V(3); % target Z-velocity in gcef (m)

    % Calculate the distance at this time, and additional expressions to
    % simplify the calculations below
    R = norm([X Y Z]); % distance (m)
    RdotV = dot([X Y Z], [VX VY VZ]);
    h = norm(cross([X Y Z], [VX VY VZ])); % specific angular momentum (m^2/s)

    % Iteratively calculate the acceleration of the chaser in the LVLH
    % reference frame, using the first-order approximation of the
    % linearised relative equations of motion
    dax = (2*mu/R^3 + h^2/R^4)*dx - 2*RdotV/R^4*h*dy + 2*h/R^2*dvy;
    day = -(mu/R^3 - h^2/R^4)*dy + 2*RdotV/R^4*h*dx - 2*h/R^2*dvx;
    daz = -mu/R^3*dz;
    
    % Compile the dydt vector for the next iteration
    dydt = [dvx dvy dvz dax day daz]'; % compiled vector 
    
end