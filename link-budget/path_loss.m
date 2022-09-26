%this function takes an input of distance (km) and frequency (Hz),
%returning the path loss (dBi) due to flux spread
function L_FS = path_loss(d)

temp = 4*pi*d^2;
L_FS = 10*log10(temp);

end
