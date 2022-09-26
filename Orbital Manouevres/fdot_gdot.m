% Lagrange coefficients derivaties calculator in terms of the change in
% universal anomaly

function [fdot, gdot] = fdot_gdot(x, r, r_0, a)

    % Global Parameters
    mu = 3.986004418e+14; % standard gravitational parameter for Earth (m^3/s^2)   

    z = a*x^2; 
    
    fdot = sqrt(mu)/r/r_0*(z*stumpffS(z) - 1)*x; 
    gdot = 1 - x^2/r*stumpffC(z);






end