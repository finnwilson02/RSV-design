function L_pol = loss_polarisation(theta)

temp = 0.5 + 0.5*cosd(2*theta);
L_pol = -10*log10(temp);

if L_pol > 41
    L_pol = 41;
end

end