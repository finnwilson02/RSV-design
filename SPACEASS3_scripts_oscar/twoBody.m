function dydt = twoBody(t, f)

    % This function calculates the velocities and accelerations of a
    % satellite for a two-body problem with the Earth. 
    mu = 3.986e+14; % standard gravitational parameter for Earth (m^3/s^2)
    
    % Extract current position and velocity
    x = f(1);
    y = f(2);
    z = f(3);
    vx = f(4);
    vy = f(5);
    vz = f(6);
    
    % Find the displacement
    r_vec = [x; y; z];
    v_vec = [vx; vy; vz];
    r = norm(r_vec);
    
    r_ddot = -mu/r^3*r_vec;

    % Compile the new dydt vector, replacing the previous displacement and
    % velocity for the current velocity and acceleration respectively
    dydt = [v_vec; r_ddot];

end