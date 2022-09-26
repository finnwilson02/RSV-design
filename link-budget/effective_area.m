%this function finds the effective antenna area as a function of frequency
function Ae = effective_area(f)

c = 297e6;
temp = 4*pi/((c/f)^2);
Ae = -10*log10(temp);

end