function [P, m_total, P_mean] = schmidt_analysis(theta, V_total, V_exp, V_comp, params)
    % Schmidt analysis for Stirling engine pressure calculation

    T_h = params.hotTemperature;
    T_c = params.coldTemperature;
    T_r = params.regeneratorTemperature;
    V_reg = params.regeneratorVolume;
    R = params.gasConstant;

    % Find BDC conditions
    [~, bdc_idx] = min(abs(theta));
    V_exp_bdc = V_exp(bdc_idx);
    V_comp_bdc = V_comp(bdc_idx);
    P_bdc = params.pressureAtBDC;

    % Calculate total mass
    denominator_bdc = V_comp_bdc/T_c + V_reg/T_r + V_exp_bdc/T_h;
    m_total = P_bdc * denominator_bdc / R;

    % Calculate pressure throughout cycle
    denominator = V_comp./T_c + V_reg/T_r + V_exp./T_h;
    P = (m_total * R) ./ denominator;

    if any(P <= 0)
        error('Negative or zero pressure detected');
    end

    P_mean = mean(P);
end