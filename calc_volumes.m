function [V_total, V_exp, V_comp, x_power, x_disp] = calc_volumes(theta, params)
    % CALC_VOLUMES Calculate instantaneous volumes for Beta-type Stirling engine
    % 
    % Inputs:
    %   theta - crank angle array (radians)
    %   params - engine parameters structure
    %
    % Outputs:
    %   V_total - total gas volume (m^3)
    %   V_exp - expansion space volume (m^3)
    %   V_comp - compression space volume (m^3)
    %   x_power - power piston position (m)
    %   x_disp - displacer position (m)
    
    % Extract parameters for clarity
    r_p = params.power.crank_length;      % Power piston crank radius
    l_p = params.power.rod_length;        % Power piston connecting rod
    r_d = params.disp.crank_length;       % Displacer crank radius
    l_d = params.disp.rod_length;         % Displacer connecting rod
    phase = params.phase_shift;           % Phase angle between pistons
    A = params.cylinder.area;             % Cylinder cross-sectional area
    
    % Validate crank-slider geometry
    if l_p <= r_p
        error('Power piston connecting rod must be longer than crank radius');
    end
    if l_d <= r_d
        error('Displacer connecting rod must be longer than crank radius');
    end
    
    % Calculate piston positions using crank-slider kinematics
    % Position measured from TDC (Top Dead Center)
    
    % Power piston position (compression space)
    beta_p = asin(r_p * sin(theta) / l_p);  % Connecting rod angle
    x_power = r_p * cos(theta) + l_p * cos(beta_p);
    x_power_tdc = r_p + l_p;  % Position at TDC
    x_power = x_power_tdc - x_power;  % Distance from TDC
    
    % Displacer position (with phase shift)
    theta_disp = theta - phase;  % Displacer leads by phase angle
    beta_d = asin(r_d * sin(theta_disp) / l_d);
    x_disp = r_d * cos(theta_disp) + l_d * cos(beta_d);
    x_disp_tdc = r_d + l_d;
    x_disp = x_disp_tdc - x_disp;
    
    % Calculate instantaneous volumes
    
    % Compression space volume (cold space)
    % Volume between power piston and bottom of displacer
    V_comp = A * x_power + params.V_dead_cold;
    
    % Expansion space volume (hot space)
    % Volume above the displacer
    V_exp_swept = A * x_disp;
    V_exp = V_exp_swept + params.V_dead_hot;
    
    % Regenerator volume (constant)
    V_reg = params.V_regenerator;
    
    % Total gas volume
    V_total = V_comp + V_exp + V_reg;
    
    % Validate results
    if any(V_comp < 0) || any(V_exp < 0)
        error('Negative volume detected - check geometry parameters');
    end
    
    if any(V_total < params.V_dead_total)
        error('Total volume less than dead volume - check calculations');
    end
    
    % Calculate and display volume characteristics (only on first call)
    persistent displayed;
    if isempty(displayed)
        V_max_calc = max(V_total);
        V_min_calc = min(V_total);
        CR_actual = V_max_calc / V_min_calc;
        
        fprintf('Volume Analysis:\n');
        fprintf('  Max Volume: %.4f L\n', V_max_calc * 1000);
        fprintf('  Min Volume: %.4f L\n', V_min_calc * 1000);
        fprintf('  Compression Ratio: %.3f\n', CR_actual);
        fprintf('  Swept Volume (Power): %.4f L\n', params.V_swept_power * 1000);
        fprintf('  Swept Volume (Displacer): %.4f L\n', params.V_swept_disp * 1000);
        fprintf('  Dead Volume: %.4f L\n', params.V_dead_total * 1000);
        fprintf('\n');
        
        displayed = true;
    end
end