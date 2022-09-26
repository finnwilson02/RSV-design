% State vector calculator, for given initial state vector and elapsed time


function [R, V] = rv_from_r0v0(R0, V0, t)

    % Inputs are in GCEF 

    % Global Parameters
    mu = 3.986004418e+14; % standard gravitational parameter for Earth (m^3/s^2)

    % Calculate the norms of the initial vectors
    r0 = norm(R0); % (m)
    v0 = norm(V0); % (m/s)
    
    % Calculate initial radial velocity
    vr0 = dot(R0, V0)/r0; % initial radial velocity (m/s)
    
    % Calculate the semi major axis' reciprocal
    alpha = 2/r0 - v0^2/mu; % reciprocal semi major axis (1/km)
    
    % Calculate the universal anomaly 
    x = U_KeplerSolve(t, r0, vr0, alpha); % universal anomaly (km^0.5)
    
    % Obtain the Lagrange coefficients
    [f, g] = fg(x, t, r0, alpha); % Lagrange coefficients 
    
    % Calculate the final position vector in gcef
    R = f*R0 + g*V0; % final position in gcef (m)
    
    % Finally, calculate the derivatives of f and g to find the final
    % velocity vector
    r = norm(R); 
    [fdot, gdot] = fdot_gdot(x, r, r0, alpha); % Lagrange coefficient first derivatives
    
    V = fdot*R0 + gdot*V0; % final velocity in gcef (m/s)
    

end
