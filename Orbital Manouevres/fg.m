% Lagrange coefficients f and g calculator in terms of the change in
% universal anomaly

function [f, g] = fg(x, t, r_0, a)

    % x = universal anomaly 
    % t = time elapsed since initial vector
    % r_0 = initial radial position 
    % a = reciprocal of semi-major axis

    % Global Parameters
    mu = 3.986004418e+14; % standard gravitational parameter for Earth (m^3/s^2)    
    
    % Calculate z to find the Stumpff functions
    z = a*x^2; 
    
    % Calculate the output Lagrange coefficients 
    f = 1 - x^2/r_0*stumpffC(z); 
    g = t - 1/sqrt(mu)*x^3*stumpffS(z);


end