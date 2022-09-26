%this function takes distance from camera to satellite r (m) and returns
%required angular resolution theta (degrees) of the camera. the use can
%input print=1 to display results to console, or print=0 to return result
%with no message
function fov = angular_resolution(r,print)

sat_radius = 8.1; %distance from center of satellite to edge (m)
required_detail = 2.5e-3;   %required detail on satellite (m)

r1 = sqrt((sat_radius-required_detail)^2 + r^2); 
r2 = sqrt(sat_radius^2 + r^2);

num = r1^2 + r2^2 - required_detail^2;
den = 2*r1*r2;
theta = acosd(num/den);

fov = 0:0.1:100; %field of view of camera (degrees)
xsize = 4096; %number of pixels in x direction
ysize = 3072; %number of pixels in y direction

fovx = theta*xsize;
fovy = theta*ysize;

fov = [fovx, fovy];

if print == 1
    
    fprintf("pixel size = [%f, %f] degrees\n",theta_x,theta_y);
    fprintf("required angular resolution for %fmm pixels: %f degrees\n",required_detail*1e3,theta);

    if theta > theta_x && theta > theta_y

        fprintf("sufficient detail is acheived at this distance\n");
    else
        fprintf("resolution is insufficient to provide required detail at this distance\n");
    end
end
end