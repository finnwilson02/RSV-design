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
% OAP connected to a 1300 kg satellite.

% Constant declaration
g0 = 9.81;          % gravity at sea level [m/s^2]

Isp = 3000;         % ion thruster specific impulse [s]
delta_v = 338;     % delta-v requirement for mission lifetime [m/s]

m_sat = 1300;       % satellite dry mass [kg]
m_oap_comp = 613.4;    % oap component mass [kg]
m_per_a_str = 10;   % mass per area of structure [kg/m^2]
den_xe = 1500;       % xenon assumed density [kg/m^3]

m_prop_prev = 0;    % previous propellant mass [kg]
m_str = 200;        % initial assumed structural mass [kg]

% Iteratively calculate the propellant mass and structural mass until the
% propellant mass converges, at which point the loop breaks.
while true
    
    % Find the total final mass of the OAP and satellite
    m_f = m_sat+m_oap_comp+m_str;
    
    % Find the propellant mass using ideal rocket equation
    m_prop = 1.1*m_f*(exp(delta_v/(Isp*g0))-1);
    
    % Determine if the propellant mass has converged to required tolerance
    if round(m_prop_prev,2) == round(m_prop,2)
        break;
    else
        m_prop_prev = m_prop;
    end
    
    % Find volume of xenon propellant required
    v_xe = m_prop/den_xe;
    
    % Add 10% volume to structure to account for components/inefficiencies.
    v_str = 1.1*v_xe;
    
    % Calculate value of L using symbolic notation.
    syms l
    eqn = l^3+4*l^2+4*l == v_str;
    l_sym = vpasolve(eqn, l);
    
    % Convert to numerical storage, taking the first value since it will be real.
    l_num = double(l_sym(1));
    a_str = 6*l_num^2+28*l_num+32;
    
    % Calculate mass of OAP from surface area
    m_str = m_per_a_str*a_str;
    
end

% Display relevant values to console
disp(m_prop);
disp(m_str);
disp(l_num);

