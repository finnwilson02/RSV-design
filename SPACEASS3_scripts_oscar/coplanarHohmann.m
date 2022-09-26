% This function calculates the required delta-v's for a coplanar Hohmann
% transfer (w/ common apse line), burning at the apogee of the first orbit
% and for increasing altitude

function [delta_vs, hohmann] = coplanarHohmann(a1, e1, a3, e3, omega1, omega3, apogeeBurn)
    
% Global parameters
    mu = 3.986004418e+14; % standard gravitational parameter for Earth (m^3/s^2)

    % Orbit 1 is initial orbit
    % Orbit 2 is the transfer orbit
    % Orbit 3 is the target orbit

    % Calculate the values needed at point A (first burn) and point B
    % (second burn)
    
    r_a1 = a1*(1 + e1); 
    r_p1 = a1*(1 - e1);
    r_a3 = a3*(1 + e3); 
    r_p3 = a3*(1 - e3);
    
    % Check if the boolean value for apogeeBurn is true of false
    % If true, the first burn is at apogee of orbit 1
    if apogeeBurn == true
        rA = r_a1; % (m)
        
        % Burn at apogee of 1, check point B
        % Check the common apse line if there is a difference in omega
        if omega1 == omega3
            rB = r_p3; % (m)
        else
            rB = r_a3; % (m)
        end
        
    else
        rA = r_p1; % (m)
        
        % Burn at perigee of 1, check point B
        % Check the common apse line if there is a difference in omega
        if omega1 == omega3
            rB = r_a3; % (m)
        else
            rB = r_p3; % (m)
        end
    end
    
    % Calculate the geometry of the Hohmann ellipse
    e2 = abs((rB - rA)/(rB + rA)); % eccentricity
    a2 = (rB + rA)/2; % semi-major axis (m)
    h2 = sqrt(mu*a2*(1 - e2^2)); % specific angular momentum of transfer orbit (m^2/s)
    
    % Calculate the required delta-v for the Hohmann at each position
    vA_2 = h2/rA; % vel for Hohmann tranfer at point A (m/s)
    vB_2 = h2/rB; % vel for Hohmann tranfer at point B (m/s)
    
    % Calculate the velocity of orbit 1 point A and orbit 3 point B
    h1 = sqrt(mu*a1*(1 - e1^2)); % specific angular momentum of transfer orbit (m^2/s)
    h3 = sqrt(mu*a3*(1 - e3^2)); % specific angular momentum of transfer orbit (m^2/s)
    
    vA_1 = h1/rA; % vel for orbit 1 at point A (m/s)
    vB_3 = h3/rB; % vel for orbit 3 at point B (m/s)

    % Calculate the delta-v's required
    delta_vA = vA_2 - vA_1; % burn A delta-v (m/s)
    delta_vB = vB_3 - vB_2; % burn B delta-v (m/s)

    % Compile the results into a row vector
    delta_vs = [delta_vA; delta_vB]; % required delta-vs
    % Compile the Hohmann elements
    hohmann = [h2, e2]; 

end