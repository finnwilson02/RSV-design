% Convert from GCEF to LVLH reference frame, given two state vectors
% This function outputs column vectors in LVLH, and takes row vectors

function [rrel_x, vrel_x, arel_x] = gcef_to_LVLH(rA, vA, rB, vB)

    % Global parameters
    mu = 3.986004418e+14; % standard gravitational parameter for Earth (m^3/s^2)

    % Body A is the target
    % Body B is the chaser
    % Inputs are in SI units (m, m/s)
    % Inputs are row vectors 
    
    % Create the basis of the LVLH reference frame centred at the target
    hA = cross(rA, vA); % angular momentum of A (m^2/s)

    i = rA/norm(rA); % X-direction unit vector for LVLH
    k = hA/norm(hA); % Z-direction unit vector for LVLH
    j = cross(k, i); % Y-direction unit vector for LVLH
    
    % Create the rotation operator for gcef to LVLH conversion
    Q_Xx = [i; j; k]; % rotation operator

    % Calculate the LVLH frame's angular velocity and acceleration vectors
    Omega = hA/(norm(rA)^2); % angular momentum of LVLH frame (rad/s)
    Omegadot = -2*((dot(vA, rA))/(norm(rA)^2))*Omega; % angular acceleration of LVLH frame (rad/s/s)
    
    % Calculate the absolute accelerations of A and B using the two-body
    % relation
    aA = -mu*rA/(norm(rA)^3); % absolute acceleration of A (m/s/s)
    aB = -mu*rB/(norm(rB)^3); % absolute acceleration of A (m/s/s)    

    % Now, calculate the position, velocity and acceleration of B relative
    % to A in gcef
    rrel = rB - rA; % relative position in gcef (m)
    vrel = vB - vA - cross(Omega, rrel); % relative velocity in gcef (m)
    arel = aB - aA - cross(Omegadot, rrel) - cross(Omega, cross(Omega, rrel)) - 2*cross(Omega, vrel); % relative acceleration in gcef (m/s/s)
    
    
    % Finally, convert the relative vectors from gcef to LVLH using the
    % rotation operator
    rrel_x = (Q_Xx*(rrel')); % LVLH position of B (m)
    vrel_x = (Q_Xx*(vrel')); % LVLH position of B (m)
    arel_x = (Q_Xx*(arel')); % LVLH position of B (m)
    
    

end