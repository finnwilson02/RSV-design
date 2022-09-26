%this function calculates the gain of a parabolic antenna, from inputs
%efficiency factor n, diameter of dish d (m), and uplink frequency f (Hz)
function G = gain_parabolic(n,d,f)

c = 3e8; %speed of light (m/s^2)
temp = n*((pi*d*f)/c)^2;
G = 10*log10(temp); %gain (dBi)

end