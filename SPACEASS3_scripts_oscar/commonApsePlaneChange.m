% This function calculates the required delta-v to be applied at the planar
% intersection of two orbits, being orbit 1, the initial orbit, and orbit
% 2, the target orbit. This calculation is ONLY for apogee planar change,
% with common apse line. 


function [delta_v] = commonApsePlaneChange(v1, v2, delta)

    % 1 is the initial orbit
    % 2 is the target
    delta_v = sqrt((v2 - v1)^2 + v1*v2*(sin(delta/2))^2); % required delta-v (m/s)
    
    



end