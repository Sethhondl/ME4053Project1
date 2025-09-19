function [P, m_total, P_mean] = schmidt_analysis(theta, V_total, V_exp, V_comp, params)
%SCHMIDT_ANALYSIS Calculate pressure using isothermal Schmidt theory
%   Assumes instantaneous pressure equilibrium and isothermal spaces

    %% ========== EXTRACT PARAMETERS ==========

    T_h = params.hotTemperature;       % K - Hot space temperature
    T_c = params.coldTemperature;      % K - Cold space temperature
    T_r = params.regeneratorTemperature;  % K - Regenerator temperature
    V_reg = params.regeneratorVolume;  % m³ - Regenerator dead volume
    R = params.gasConstant;            % J/(kg·K) - Gas constant for air


    %% ========== FIND REFERENCE CONDITIONS ==========

    % Locate bottom dead center (BDC) position - maximum volume
    [~, bdc_idx] = max(V_total);

    % Get volumes at BDC
    V_exp_bdc = V_exp(bdc_idx);
    V_comp_bdc = V_comp(bdc_idx);
    P_bdc = params.pressureAtBDC;


    %% ========== CALCULATE TOTAL GAS MASS ==========

    % Apply ideal gas law at BDC condition
    denominator_bdc = V_comp_bdc/T_c + V_reg/T_r + V_exp_bdc/T_h;
    m_total = P_bdc * denominator_bdc / R;  % kg


    %% ========== CALCULATE INSTANTANEOUS PRESSURE ==========

    % Schmidt equation: P = m_total * R / (V_c/T_c + V_r/T_r + V_h/T_h)
    denominator = V_comp./T_c + V_reg/T_r + V_exp./T_h;
    P = (m_total * R) ./ denominator;  % Pa


    %% ========== VALIDATION ==========

    if any(P <= 0)
        error('Negative or zero pressure detected - check volume calculations');
    end

    % Calculate mean pressure for reporting
    P_mean = mean(P);

end