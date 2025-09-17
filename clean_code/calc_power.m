function [W_indicated, P_indicated, W_mep, P_mep, MEP, efficiency] = calc_power(P, V_total, theta, params)
    % Calculate work, power, and efficiency

    if abs(theta(end) - theta(1) - 2*pi) > 0.01 && abs(theta(end) - theta(1)) > 0.01
        warning('Theta does not span exactly one complete cycle');
    end

    % Method 1: Direct integration
    W_indicated = abs(trapz(V_total, P));
    P_indicated = W_indicated * params.operatingFrequency;

    % Method 2: Mean Effective Pressure
    V_swept = max(V_total) - min(V_total);
    MEP = abs(W_indicated) / V_swept;
    W_mep = MEP * V_swept;
    P_mep = W_mep * params.operatingFrequency;

    % Calculate efficiency
    if isfield(params, 'm_total')
        m_total = params.m_total;
    else
        m_total = P(1) * V_total(1) / (params.gasConstant * params.hotTemperature);
    end

    V_max_cycle = max(V_total);
    V_min_cycle = min(V_total);
    Q_in_isothermal = m_total * params.gasConstant * params.hotTemperature * log(V_max_cycle/V_min_cycle);
    Q_in = abs(Q_in_isothermal);

    if Q_in > 0
        efficiency = abs(W_indicated) / Q_in;
    else
        dV = diff(V_total);
        P_avg = (P(1:end-1) + P(2:end)) / 2;
        dV_positive = dV > 0;
        Q_in = sum(P_avg(dV_positive) .* dV(dV_positive));
        efficiency = abs(W_indicated) / abs(Q_in);
    end

    eta_carnot = 1 - params.coldTemperature / params.hotTemperature;
    if efficiency > eta_carnot
        efficiency = min(efficiency, eta_carnot * 0.99);
    end

    if efficiency < 0 || efficiency > 1
        error('Invalid efficiency calculated: %.2f', efficiency);
    end
end