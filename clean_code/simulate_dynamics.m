function [omega, alpha, rpm, Cs_actual] = simulate_dynamics(T_total, theta, I_flywheel, params)
    omega_avg = params.averageRPM * 2 * pi / 60;
    T_mean = mean(T_total);

    omega = zeros(size(theta));
    alpha = zeros(size(theta));
    omega(1) = omega_avg;

    for i = 2:length(theta)
        dtheta = theta(i) - theta(i-1);
        dt = dtheta / omega(i-1);

        T_net = T_total(i-1) - T_mean;
        alpha(i-1) = T_net / I_flywheel;
        omega(i) = omega(i-1) + alpha(i-1) * dt;

        if omega(i) <= 0
            error('Angular velocity became negative - increase flywheel inertia');
        end
    end

    alpha(end) = alpha(end-1);
    rpm = omega * 60 / (2*pi);

    omega_max = max(omega);
    omega_min = min(omega);
    omega_mean = mean(omega);
    Cs_actual = (omega_max - omega_min) / omega_mean;

    if Cs_actual > params.flywheelCoefficientOfFluctuation * 1.1
        warning('Actual Cs (%.4f) exceeds target (%.4f)', Cs_actual, params.flywheelCoefficientOfFluctuation);
    end
end