% Orbital elements to state vector converter

function [r, v] = coe2stateVecs(coe)

    % coe is the orbital elements, inputted as [h; e; Omega; i; omega; theta]
    % Inputs in SI and radians
    
    % Global parameters
    mu = 3.986004418e+14; % standard gravitational parameter for Earth (m^3/s^2)

    % Extract classical orbital elements
    h = coe(1); % specific angular momentum (m^2/s)
    e = coe(2); % eccentricity 
    Omega = coe(3); % right asc of ascending node (rad)
    i = coe(4); % inclination (rad)
    omega = coe(5); % arg of perigee (rad)
    theta = coe(6); % true anomaly (rad)
    
    % Convert the orbital elements to perifocal frame state vectors 
    r_p = (h^2/mu)*(1/(1+ e*cos(theta)))*[cos(theta); sin(theta); 0]; % position (m)
    v_p = (mu/h)*[-sin(theta); e + cos(theta); 0]; % velocity (m)
    
    % Construct the rotation operator
    O = Omega;
    w = omega;
    Q_xX = [-sin(O)*cos(i)*sin(w) + cos(O)*cos(w), -sin(O)*cos(i)*cos(w) - cos(O)*sin(w), sin(O)*sin(i);
        cos(O)*cos(i)*sin(w) + sin(O)*cos(w), cos(O)*cos(i)*cos(w) - sin(O)*sin(w), -cos(O)*sin(i);
        sin(i)*sin(w), sin(i)*cos(w), cos(i)]; % perifocal to GCEF rotation operator
    
    % Apply the transformation to each vector 
    r = (Q_xX*r_p)'; % gcef position (m)
    v = (Q_xX*v_p)'; % gcef velocity (m)
    
end
    

