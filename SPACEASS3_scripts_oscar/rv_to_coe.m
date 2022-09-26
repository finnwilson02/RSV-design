%

function elements = rv_to_coe(R, V)

    % Define the standard gravitational parameter
    mu = 3.986004418e+14; % standard gravitational parameter for Earth (m^3/s^2)

    % Find the magnitude of each state vector
    r = norm(R);
    v = norm(V);
    
    % Calculate the radial velocity
    vr = dot(R, V)/r; % radial velocity (m/s)
    
    % Calculate specific angular momentum
    H = cross(R ,V); % m^2/s
    h = norm(H); % m^2/s
    
    % Calculate the inclination angle
    i = acos(H(3)/h); % inclination (rad)
    
    % Calculate nodal line vector 
    N = cross([0; 0; 1], H);
    n = norm(N);
    
    % Calculate the ascending node's right ascension
    if n ~= 0 
        
        Omega = acos(N(1)/n); % right ascension of asc node (rad)
        
        if N(2) < 0
            
            Omega = 2*pi - Omega; % right ascension of asc node, expressed within [0, 2pi] (rad)
        end
    else
        Omega = 0;
        
    end
    
    % Find the eccentricity vector to calculate the remainding elements
    E = (1/mu)*((v^2 - mu/r)*R - r*vr*V); % eccentricity vector
    e = norm(E); % magnitude of E
    
    % Calculate the argument of perigee with this information, accounting
    % for quadrant ambiguity
    if n ~= 0 
        
        if e > eps
            
            omega = acos(dot(N, E)/n/e); % arg of perigee (rad)
    
            if E(3) < 0
                
                omega = 2*pi - omega; % reducing to the correct range (rad)
                
            end
        else
            omega = 0; % if n is 0, the argument of perigee lies on the equatorial plane
        end
    else
        
        omega = 0; % rad
        
    end
    
    % Now, calculate the true anomaly accounting for quadrant ambiguity
    if e > eps
        
        theta = real(acos(dot(E, R)/e/r)); % true anomaly (rad)
        
        if vr < 0
            
            theta = 2*pi - theta; % checking and converting to within 0, 2pi
            
        end
    else
        % check the cross product of the nodal line and position vector
        C = cross(N, R); 
        
        % if the third term is positive or negative, calculcate the true
        % anomaly for each
        if C(3) >= 0
            
            theta = acos(dot(N, R)/n/r); % true anomaly (rad)
        else
            theta = 2*pi - acos(dot(N, R)/n/r); % (rad)
        end
        
        
    end
    
    % Finally, calculate the semi-major axis
    a = h^2/mu/(1 - e^2); % semi-major axis (m)
    
    % Compile the output vector full of the elements
    elements = [a, e, h, i, theta, Omega, omega]; 
    





end