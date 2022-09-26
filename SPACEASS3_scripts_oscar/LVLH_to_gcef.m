% This function converts from the LVLH rotating frame to the GCEF reference
% frame

function [R, V] = LVLH_to_gcef(COE, rrel_x, vrel_x)

    % Global parameters
    mu = 3.986004418e+14; % standard gravitational parameter for Earth (m^3/s^2)

    % COE = classical orbital elements for the target orbit, A
    
    % Firstly, extract the state vector of the target in GCEF frame
    [rA_X, vA_X] = coe2stateVecs(COE);
    
    % Create the basis of the LVLH reference frame centred at the target
    hA = cross(rA_X, vA_X); % angular momentum of A (m^2/s)

    i = rA_X/norm(rA_X); % X-direction unit vector for LVLH
    k = hA/norm(hA); % Z-direction unit vector for LVLH
    j = cross(k, i); % Y-direction unit vector for LVLH
    
    % Create the rotation operator for LVLH to GCEF conversion
    Q_xX = [i; j; k]'; % rotation operator from LVLH to GCEF

    % Calculate the LVLH frame's angular velocity and acceleration vectors
    Omega = hA/(norm(rA_X)^2); % angular momentum of LVLH frame (rad/s)
    
    % Convert from LVLH relative motion to GCEF relative motion of
    % spacecraft B
    rrel_X = (Q_xX*(rrel_x)); % GCEF position of B relative to A (m)
    vrel_X = (Q_xX*(vrel_x)); % GCEF velocity of B relative to A (m)
    
    % Now, the GCEF relative motion vectors must be converted to absolute
    % motion 
    R = rrel_X' + rA_X; % the absolute position of the chaser, B (m)
    V = vrel_X' + vA_X + (cross(Omega,rrel_X)); % the absolute velocity of the chaser, B (m/s)
    
end