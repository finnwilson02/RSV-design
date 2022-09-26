function dydt = linearRelMotionRatesTEST(t, y, R0, V0)

    % Global parameters
    mu = 3.986004418e+14; % standard gravitational parameter for Earth (m^3/s^2)
    J2 = 0.00108263; % second zonal harmonics effect 
    R_e = 6.3781e+6; % Mean radius of Earth (m)
    
    % Extract the input vector's variables
    dx = y(1); % x-pos of chaser (m)
    dy = y(2); % y-pos of chaser (m)
    dz = y(3); % z-pos of chaser (m)
    dvx = y(4); % x-vel of chaser (m)
    dvy = y(5); % x-vel of chaser (m)
    dvz = y(6); % x-vel of chaser (m)
    
    % Calculate the J2's effect on the orbit
    X = y(7);
    Y = y(8);
    Z = y(9);
    VX = y(10); % target X-velocity in gcef (m)
    VY = y(11); % target Y-velocity in gcef (m)
    VZ = y(12); % target Z-velocity in gcef (m)    
    
    % Find the displacement
    R = norm([X Y Z]);
    
    % Calculate the acceleration at the current time
    ax = -mu*X/R^3 + (3/2)*((J2*mu*R_e^2)/(R^4))*(X/R)*((5*Z^2)/(R^2) - 1); % (m/s/s)
    ay = -mu*X/R^3 + (3/2)*((J2*mu*R_e^2)/(R^4))*(X/R)*((5*Z^2)/(R^2) - 1); % (m/s/s)
    az = -mu*X/R^3 + (3/2)*((J2*mu*R_e^2)/(R^4))*(X/R)*((5*Z^2)/(R^2) - 3); % (m/s/s)
    
    % Calculate the distance at this time, and additional expressions to
    % simplify the calculations below
    RdotV = dot([X Y Z], [VX VY VZ]);
    h = norm(cross([X Y Z], [VX VY VZ])); % specific angular momentum (m^2/s)

    % Iteratively calculate the acceleration of the chaser in the LVLH
    % reference frame, using the first-order approximation of the
    % linearised relative equations of motion
    dax = (2*mu/R^3 + h^2/R^4)*dx - 2*RdotV/R^4*h*dy + 2*h/R^2*dvy;
    day = -(mu/R^3 - h^2/R^4)*dy + 2*RdotV/R^4*h*dx - 2*h/R^2*dvx;
    daz = -mu/R^3*dz;
    
    % Compile the dydt vector for the next iteration
    dydt = [dvx dvy dvz dax day daz VX VY VZ ax ay az]'; % compiled vector 
    
end