

function dum = Phi_vr(n, t)

    dum = [3*n*sin(n*t), 0, 0; 
    6*n*(cos(n*t) - 1), 0, 0; 
    0, 0, -n*sin(n*t)];

end