function [P, m_total, P_mean] = schmidt_analysis(theta, V_total, V_exp, V_comp, params)
    % SCHMIDT_ANALYSIS Calculate pressure using Schmidt theory for Stirling engine
    %
    % Inputs:
    %   theta - crank angle array (radians)
    %   V_total - total gas volume array (m^3)
    %   V_exp - expansion space volume array (m^3)
    %   V_comp - compression space volume array (m^3)
    %   params - engine parameters structure
    %
    % Outputs:
    %   P - instantaneous pressure array (Pa)
    %   m_total - total mass of working fluid (kg)
    %   P_mean - mean cycle pressure (Pa)
    
    % Extract temperatures
    T_h = params.T_hot;         % Hot space temperature
    T_c = params.T_cold;        % Cold space temperature
    T_r = params.T_regenerator; % Regenerator temperature
    
    % Extract volumes
    V_reg = params.V_regenerator;  % Regenerator volume (constant)
    
    % Gas constant
    R = params.gas.R;  % Specific gas constant for working fluid
    
    % Schmidt analysis assumes isothermal spaces and ideal regenerator
    % The pressure is uniform throughout the engine at any instant
    
    % Calculate mass distribution (isothermal assumption)
    % P*V = m*R*T, so m = P*V/(R*T)
    % Total mass is conserved: m_total = m_hot + m_cold + m_reg
    
    % For Schmidt analysis, pressure is given by:
    % P = (m_total * R) / (V_comp/T_c + V_reg/T_r + V_exp/T_h)
    
    % First, find total mass from initial conditions (at BDC)
    % Assuming theta=0 is at BDC for power piston
    [~, bdc_idx] = min(abs(theta));
    V_total_bdc = V_total(bdc_idx);
    V_exp_bdc = V_exp(bdc_idx);
    V_comp_bdc = V_comp(bdc_idx);
    P_bdc = params.P_bdc;
    
    % Calculate total mass from BDC conditions
    denominator_bdc = V_comp_bdc/T_c + V_reg/T_r + V_exp_bdc/T_h;
    m_total = P_bdc * denominator_bdc / R;

    % Now calculate pressure throughout the cycle
    denominator = V_comp./T_c + V_reg/T_r + V_exp./T_h;
    P = (m_total * R) ./ denominator;

    % Validate mass conservation (engineering assumption)
    % Calculate mass in each space throughout cycle
    m_comp = P .* V_comp / (R * T_c);
    m_exp = P .* V_exp / (R * T_h);
    m_reg = P .* V_reg / (R * T_r);
    m_cycle_total = m_comp + m_exp + m_reg;

    % Check mass conservation
    mass_variation = (max(m_cycle_total) - min(m_cycle_total)) / m_total;
    if mass_variation > 1e-6  % Allow for numerical precision
        warning('Mass conservation violation detected: %.2e relative variation', mass_variation);
        fprintf('  Total mass should be constant at %.6f g\n', m_total * 1000);
        fprintf('  But varies between %.6f and %.6f g\n', ...
                min(m_cycle_total) * 1000, max(m_cycle_total) * 1000);
    end
    
    % Validate pressure
    if any(P <= 0)
        error('Negative or zero pressure detected');
    end
    
    if any(P > 100e6)  % 100 MPa limit
        warning('Extremely high pressure detected (>100 MPa)');
    end
    
    % Calculate mean pressure
    P_mean = mean(P);
    
    % Calculate pressure characteristics
    P_max = max(P);
    P_min = min(P);
    pressure_ratio = P_max / P_min;
    
    % Display pressure analysis (only on first call)
    persistent displayed;
    if isempty(displayed)
        fprintf('Pressure Analysis (Schmidt Theory):\n');
        fprintf('  Total Working Fluid Mass: %.6f g\n', m_total * 1000);
        fprintf('  Maximum Pressure: %.3f MPa\n', P_max / 1e6);
        fprintf('  Minimum Pressure: %.3f MPa\n', P_min / 1e6);
        fprintf('  Mean Pressure: %.3f MPa\n', P_mean / 1e6);
        fprintf('  Pressure Ratio: %.3f\n', pressure_ratio);
        fprintf('\n');
        
        displayed = true;
    end
    
    % Additional validation
    if pressure_ratio > 10
        warning('Very high pressure ratio (>10) - may indicate parameter issues');
    end
end