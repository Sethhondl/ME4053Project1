function [P, m_total, P_mean] = schmidt_analysis(theta, V_total, V_exp, V_comp, params)
    T_h = params.hotTemperature;
    T_c = params.coldTemperature;
    T_r = params.regeneratorTemperature;
    V_reg = params.regeneratorVolume;
    R = params.gasConstant;

    [~, bdc_idx] = max(V_total);
    V_total_bdc = V_total(bdc_idx);
    V_exp_bdc = V_exp(bdc_idx);
    V_comp_bdc = V_comp(bdc_idx);
    P_bdc = params.pressureAtBDC;

    denominator_bdc = V_comp_bdc/T_c + V_reg/T_r + V_exp_bdc/T_h;
    m_total = P_bdc * denominator_bdc / R;

    denominator = V_comp./T_c + V_reg/T_r + V_exp./T_h;
    P = (m_total * R) ./ denominator;

    P_mean = mean(P);

    if any(P < 0)
        error('Negative pressure detected');
    end

    if max(P) > 100e6
        warning('Extremely high pressure detected');
    end
end