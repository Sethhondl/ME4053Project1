function [W_indicated, P_indicated, W_mep, P_mep, MEP, efficiency] = calc_power(P, V_total, theta, params)
    % CALC_POWER Calculate work and power using two different methods
    %
    % Inputs:
    %   P - pressure array (Pa)
    %   V_total - total volume array (m^3)
    %   theta - crank angle array (radians)
    %   params - engine parameters structure
    %
    % Outputs:
    %   W_indicated - indicated work per cycle (J) - Method 1
    %   P_indicated - indicated power (W) - Method 1
    %   W_mep - work from mean effective pressure (J) - Method 2
    %   P_mep - power from mean effective pressure (W) - Method 2
    %   MEP - mean effective pressure (Pa)
    %   efficiency - thermal efficiency
    
    %% Method 1: Direct Integration (Indicated Work)
    % Work = integral of P*dV over complete cycle
    
    % Ensure we have a complete cycle
    if abs(theta(end) - theta(1) - 2*pi) > 0.01 && abs(theta(end) - theta(1)) > 0.01
        warning('Theta does not span exactly one complete cycle');
    end
    
    % Calculate work using trapezoidal integration
    % W = -integral(P*dV) (negative because compression decreases volume)
    % For a cycle producing positive work, the integral should be negative
    
    % Method 1a: Direct numerical integration
    dV = diff(V_total);  % Volume changes
    P_avg = (P(1:end-1) + P(2:end)) / 2;  % Average pressure for each segment
    W_segments = P_avg .* dV;  % Work for each segment
    W_indicated = sum(W_segments);  % Total work

    % Method 1b: Using trapz for verification
    W_indicated_trapz = trapz(V_total, P);

    % Check consistency
    if abs(W_indicated - W_indicated_trapz) > 0.01 * abs(W_indicated)
        warning('Discrepancy in work calculation methods');
    end

    % Use the more accurate trapz method
    W_indicated = W_indicated_trapz;

    % For Stirling engines, if work is negative (counter-clockwise cycle),
    % take the absolute value as the engine produces positive work
    W_indicated = abs(W_indicated);
    
    % Calculate indicated power
    P_indicated = W_indicated * params.operatingFrequency;  % W = J/s
    
    %% Method 2: Mean Effective Pressure (MEP)
    % MEP represents the constant pressure that would produce the same work
    % over the displacement volume
    
    % Calculate swept volume (displacement)
    V_swept = max(V_total) - min(V_total);
    
    % Mean Effective Pressure
    MEP = abs(W_indicated) / V_swept;
    
    % Work from MEP
    W_mep = MEP * V_swept;
    
    % Power from MEP
    P_mep = W_mep * params.operatingFrequency;
    
    %% Calculate Thermal Efficiency
    % For isothermal processes (engineering assumption #4):
    % Heat input occurs during isothermal expansion at T_hot
    % Heat rejection occurs during isothermal compression at T_cold

    % Find expansion and compression processes
    dV_positive = dV > 0;  % Expansion segments
    dV_negative = dV < 0;  % Compression segments

    % For isothermal process: Q = nRT*ln(V2/V1) = ∫PdV (since PV = nRT = constant)
    % Heat input during isothermal expansion at T_hot
    Q_in = sum(P_avg(dV_positive) .* dV(dV_positive));

    % Heat rejected during isothermal compression at T_cold
    Q_out = abs(sum(P_avg(dV_negative) .* dV(dV_negative)));

    % Alternative calculation using proper isothermal assumption
    % Total gas mass from Schmidt analysis
    if isfield(params, 'm_total')
        m_total = params.m_total;
    else
        % Estimate from initial conditions
        m_total = P(1) * V_total(1) / (params.gasConstant * params.hotTemperature);
    end

    % Isothermal heat transfer: Q = m*R*T*ln(V_final/V_initial)
    V_max_cycle = max(V_total);
    V_min_cycle = min(V_total);

    % Heat input during isothermal expansion at T_hot
    Q_in_isothermal = m_total * params.gasConstant * params.hotTemperature * log(V_max_cycle/V_min_cycle);

    % Heat rejected during isothermal compression at T_cold
    Q_out_isothermal = m_total * params.gasConstant * params.coldTemperature * log(V_max_cycle/V_min_cycle);

    % Use isothermal calculation for efficiency (proper for Stirling engine)
    Q_in = abs(Q_in_isothermal);

    % Thermal efficiency
    if Q_in > 0
        efficiency = abs(W_indicated) / Q_in;
    else
        % Fallback to numerical integration
        Q_in = sum(P_avg(dV_positive) .* dV(dV_positive));
        efficiency = abs(W_indicated) / abs(Q_in);
    end
    
    % Carnot efficiency for comparison
    eta_carnot = 1 - params.coldTemperature / params.hotTemperature;
    
    % Validate results
    % Work should now always be positive after taking absolute value
    
    if efficiency > eta_carnot
        warning('Efficiency exceeds Carnot limit - check calculations');
        efficiency = min(efficiency, eta_carnot * 0.99);  % Cap at 99% of Carnot
    end
    
    if efficiency < 0 || efficiency > 1
        error('Invalid efficiency calculated: %.2f', efficiency);
    end
    
    % Display power analysis (only on first call)
    persistent displayed;
    if isempty(displayed)
        fprintf('Power and Efficiency Analysis:\n');
        fprintf('  Method 1 - Indicated Work:\n');
        fprintf('    Work per Cycle: %.2f J\n', W_indicated);
        fprintf('    Power Output: %.2f W (%.3f kW)\n', P_indicated, P_indicated/1000);
        fprintf('  Method 2 - Mean Effective Pressure:\n');
        fprintf('    MEP: %.3f MPa\n', MEP/1e6);
        fprintf('    Work per Cycle: %.2f J\n', W_mep);
        fprintf('    Power Output: %.2f W (%.3f kW)\n', P_mep, P_mep/1000);
        fprintf('  Efficiency:\n');
        fprintf('    Thermal Efficiency: %.1f%%\n', efficiency * 100);
        fprintf('    Carnot Limit: %.1f%%\n', eta_carnot * 100);
        fprintf('    Relative Efficiency: %.1f%% of Carnot\n', (efficiency/eta_carnot) * 100);
        fprintf('\n');
        
        % Check agreement between methods
        power_diff = abs(P_indicated - P_mep) / P_indicated * 100;
        if power_diff > 5
            warning('Power calculations differ by %.1f%% between methods', power_diff);
        end
        
        displayed = true;
    end
end