%this function finds antenna loss due to rain, taking inputs of 
function L_rain = rain_loss_down(theta,tau)

h_0 = 3; %0C isotherm height at 0N, 152E (km)
h_ant = 0; %antenna height, km
D_rain = (h_0 - h_ant +0.36)/sind(theta); %distance through rain (km)

%find values appropriate for given frequency (14 GHz for uplink) at 
%https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.838-3-200503-I!!PDF-E.pdf
k_h = 0.02386;
k_v = 0.02455; 
a_h = 1.1825;
a_v = 1.1216;

%rain rate for longitude and latitude of satellite (0N, 152E) taken from
%https://www.itu.int/dms_pubrec/itu-r/rec/p/R-REC-P.837-1-199408-S!!PDF-E.pdf
R = 22; %rain rate (mm/h). chosen for 99.9% availability

%calculate k and alpha (coefficients)
k = (k_h + k_v + (k_h - k_v)*cosd(theta)^2*cosd(2*tau));
a = (k_h*a_h - k_v*a_v + (k_h*a_h - k_v*a_v)*cosd(theta)^2*cosd(2*tau));
a = 1.1;

gamma = k*R^a; %attenuation/km
L_rain = gamma*D_rain; %loss due to rain 

end