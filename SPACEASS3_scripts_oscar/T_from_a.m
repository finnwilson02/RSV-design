

function T = T_from_a(a)

    % Global Parameters
    mu = 3.986004418e+14; % standard gravitational parameter for Earth (m^3/s^2)   
    
    T = 2*pi*sqrt((a^3)/mu); % period of SSO (s)

end