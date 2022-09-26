% Universal Kepler Equation Calculator

function x = U_KeplerSolve(delta_t, r_o, vr_o, a)

    % Global Parameters
    mu = 3.986004418e+14; % standard gravitational parameter for Earth (m^3/s^2)

    % Define an iteration tolerance 
    error = 1*10^-8; 
    nMax = 1000; % maximum iteration count
    
    % Calculate the initial universal anomaly
    x = sqrt(mu)*abs(a)*delta_t; 

    % Iterate with Newton's method of approximation until the solution 
    % converges within the specified margin of error
    n = 0; % initialise iteration counter
    ratio = 1; % initial ratio 
    
    while (abs(ratio) > error) && (n <= nMax)
        
        n = n + 1; % iterate the step counter
        
        % Evaluate the Stumpff functions
        z = a*x^2; % calculate z to obtain the Stumpff functions
        C = stumpffC(z); 
        S = stumpffS(z); 
        
        % Evaluate the function and its derivative in order to evaluate the
        % current iteration's x-value
        F = r_o*vr_o/sqrt(mu)*x^2*C + (1 - a*r_o)*x^3*S + r_o*x - sqrt(mu)*delta_t; % f(x)
        dFdx = r_o*vr_o/sqrt(mu)*x*(1 - a*x^2*S) + (1 - a*r_o)*x^2*C + r_o; % f'(x)
        
        % Calculate the ratio
        ratio = F/dFdx;

        % Calculate x for this iteration
        x = x - ratio; 

    end

end