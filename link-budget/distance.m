%this function finds the distance from the ground station to the satellite,
%returning d in meters from inputs of altitude (m) and angle from antenna
%to satellite (degrees)
function d = distance(r,theta)

Re = 6378100;
d = (r+Re)/sind(90 - theta);

end