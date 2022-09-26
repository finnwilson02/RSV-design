% Author details:
%   Name: Taj Wedutenko
%   SID: 500457324

% Start from a blank workspace.
clear;
clc;

% Clear existing figures and reset plotting properties
clf;
clf reset;
close all;

% This script is used to estimate the propellant and structural mass of an
% RSV.

% Constant declaration
g0 = 9.81;          % gravity at sea level [m/s^2]

ion_isp = 3000;     % ion thruster specific impulse [s]
ion_delta_v = 289.36; % ion thruster delta-v requirement for mission lifetime [m/s]

lae_isp = 325;
lae_delta_v = 2846.48284;

m_rsv_comp = 1042;    % rsv component mass [kg]
m_per_a_str = 10;   % mass per area of structure [kg/m^2]
den_ion = 1500;       % xenon assumed density [kg/m^3]
den_lae = 1370;

m_ion_prev = 0; % previous propellant mass [kg]
m_lae_prev = 0; % previous propellant mass [kg]
m_str = 1000;        % initial assumed structural mass [kg]

% Iteratively calculate the propellant mass and structural mass until the
% propellant mass converges, at which point the loop breaks.
while true
    
    % Find the total final mass of the OAP and satellite
    m_f = m_rsv_comp+m_str;
    
    % Find the ion propellant mass using ideal rocket equation
    m_ion = 1.1*m_f*(exp(ion_delta_v/(ion_isp*g0))-1);
    
    % Find the lae propellant mass using ideal rocket equation
    m_lae = 1.1*m_f*(exp(lae_delta_v/(lae_isp*g0))-1);
    
    % Determine if the propellant mass has converged to required tolerance
    if round(m_ion_prev,1) == round(m_ion,1) && ...
            round(m_lae_prev,1) == round(m_lae,1)
        break;
    else
        m_ion_prev = m_ion;
        m_lae_prev = m_lae;
    end
    
    % Find volume of ion propellant required
    v_ion = m_ion/den_ion;
    v_lae = m_lae/den_lae; 
    
    % Add 10% volume to structure to account for components/inefficiencies.
    v_str = 1.1*(v_ion+v_lae);
    
    % Calculate the value of L.
    l = v_str/32.5;
    a_str = 23*l+157;
    
    % Calculate mass of OAP from surface area
    m_str = m_per_a_str*a_str;
    
end

% Display relevant values to console
disp(m_ion);
disp(m_lae);
disp(m_str);
disp(l);

