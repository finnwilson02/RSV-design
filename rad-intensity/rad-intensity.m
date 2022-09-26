clear all;
clc;
clf;

%Parameters and Constants
deg = pi/180;   %convert degrees to radians
hours = 3600;   %convert hours to seconds
Re = 6371.1;    %radius of Earth (km)

G = 6.67408e-11;
M = 1.989e30;
mu = G*M;  %gravitational parameter (km^3 s^-2)

%Classical Orbital Elements
h = 4.45e15;   %angular momentum (kg^2/s)
e = 0.0167;  %eccentricity
RA = 360-11.26064;  %right ascension of the ascending node (deg)
incl = 0.00005; %inclination of the orbit (deg)
w = 102.94719;    %argument of perigee (deg)

Me = 100.46435;
E = kepler_E(e,Me);
TA = acosd((e - cosd(E))/(e*cosd(E)-1));  %true anomaly (deg)

%Get State Vectors
coe = [h, e, RA, incl, w, TA];  %concatenate classical orbital elements
[r0 v0] = coe2sv(coe,mu);                       %call coe2sv to get state vectors expressing orbit

%Set up timespan and call ODE function
t0 = 0;
tf = 60*60*24*365;
y0 = [r0 v0]';
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t,y] = ode45(@rates,[t0 tf],y0,options);


d = zeros(length(y),1);
for i = 1:length(y)
    
    d(i,1) = sqrt(y(i,1)^2 + y(i,2)^2 + y(i,3)^2);
end

figure(1);
plot(t/(60*60*24),d);

Hsurface = 6.33e7;
Rsun = 695e6;
I = Hsurface*(Rsun^2/d.^2);

figure(2);
plot(t/(60*60*24),I);
