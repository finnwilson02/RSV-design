% This function is made for the hohmannCoplanar calculator to plot the
% hohmann orbit, and converts the outputs to coe's which can then be
% propagated

function coeHohmann = hohmannExtract(h, e)

    % Global parameters
    mu = 3.986004418e+14; % standard gravitational parameter for Earth (m^3/s^2)

    a = h^2/(mu*(1 - e^2)); % semi-major axis (m)
    
    




end