% OAP fuel function



function propReq = propOAP(Isp, mA_0, mE, mP, dv_net, nBurns, years)

    % Isp = specific impulse of propellant
    % mA_0 = target satellite's dry mass (kg)
    % mE = OAP structural mass (kg)
    % mP = Initial estimate of optimal OAP propellant mass (kg)
    % nBurns = Amount of burns each year
    % years = years to be serviced 
    
    % Global parameters
    g0 = 9.805; % (m/s/s)
    c = g0*Isp; % exit velocity of propellant (m/s)
    
    % Define a tolerance, to tell the script whether to return from the
    % function or continue
    tol = 1; % tolerance (kg)
    
    m0 = mA_0 + mE + mP; % initial mass of the mated vehicle (kg)
    
    delta_v = dv_net/nBurns; % the required delta-v per burn (m/s)
    disp(delta_v);
    % Loop through the burns in a year to find the change in mass over the
    % year
    propReq = 0; % initialise counter for required fuel (kg)
    propReq_prev = mP; % the initial required fuel (kg)
    
    while propReq <= propReq_prev
        
        % If the fuel used in the previous iteration exceed physical
        % capabilities, return the previous iteration's result
        if m0 < (mA_0 + mE)
            propReq = propReq_prev;
            return
        end
        
        % Set the previous iteration's fuel mass to the 'previous' fuel
        % mass so it can be checked later
        propReq_prev = propReq;
        
        % Evaluate the initial mass of the spacecraft
        m0 = mA_0 + mE + propReq_prev; % (kg)
        
        % Set the previous fuel given to the spacecraft
        propReq = 0; % zero the fuel mass (kg)
        
        
        % Loop to calculate the mass after the burn durations
        for j = 1:years
            for k = 1:nBurns
                mf = m0/exp(delta_v/c); % calculate the mass after the burn (kg)
                propReq = propReq + (m0 - mf); % add the burnt fuel to the running total (kg)
                disp(propReq);
                m0 = mf; % change the mass (kg)
            end
        end
        
        
        
    end
 
        
    
end