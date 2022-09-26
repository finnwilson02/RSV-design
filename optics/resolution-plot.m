clear all;
clc;
clf;

fov = zeros(1000,2);
fovy = zeros(1000,1);

for i = 1:1000
    
    fov(i,:) = angular_resolution(i,0);
end

plot(1:1000,fov(:,1));
hold on;
plot(1:1000,fov(:,2));
axis([50,100,0,13]);
title("maximum field of view vs. distance");
xlabel("distance (m)");
ylabel("max field of view (degrees)");
legend('x FoV','y FoV');