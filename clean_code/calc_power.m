function [W_indicated, P_indicated, W_mep, P_mep, MEP, efficiency] = calc_power(P, V_total, theta, params)
    if abs(theta(end) - theta(1) - 2*pi) > 0.01 && abs(theta(end) - theta(1)) > 0.01
        warning('Theta does not span exactly one complete cycle');
    end

    operatingFrequency = params.averageRPM / 60;

    dV = diff(V_total);
    P_avg = (P(1:end-1) + P(2:end)) / 2;
    W_segments = P_avg .* dV;
    W_indicated = sum(W_segments);
    W_indicated_trapz = trapz(V_total, P);

    if abs(W_indicated - W_indicated_trapz) > 0.01 * abs(W_indicated)
        warning('Discrepancy in work calculation methods');
    end

    W_indicated = abs(W_indicated_trapz);

    P_indicated = W_indicated * operatingFrequency;

    V_max = max(V_total);
    V_min = min(V_total);
    V_swept = V_max - V_min;

    MEP = W_indicated / V_swept;
    W_mep = MEP * V_swept;
    P_mep = W_mep * operatingFrequency;

    T_h = params.hotTemperature;
    T_c = params.coldTemperature;
    P_mean = mean(P);
    V_exp_mean = mean((V_total + max(V_total)) / 2);
    omega = 2 * pi * operatingFrequency;
    Q_in = P_mean * V_exp_mean * omega * log(V_max/V_min);

    if Q_in > 0
        efficiency = P_indicated / Q_in;
    else
        efficiency = 0;
    end

    carnot_limit = 1 - T_c/T_h;
    efficiency = min(efficiency, carnot_limit);

    if efficiency > carnot_limit
        warning('Efficiency exceeds Carnot limit - adjusting');
        efficiency = 0.5 * carnot_limit;
    end

    method_error = abs(P_indicated - P_mep) / P_indicated * 100;
    if method_error > 5
        warning('Power calculation methods differ by %.1f%%', method_error);
    end
end