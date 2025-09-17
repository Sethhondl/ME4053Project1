function [omega, alpha, rpm, Cs_actual] = simulate_dynamics(T_total, theta, I_flywheel, params)
    % Simulate angular velocity variation with flywheel

    omega_avg = params.averageAngularVelocity;
    T_mean = mean(T_total);
    T_load = T_mean;
    T_net = T_total - T_load;
    alpha = T_net / I_flywheel;

    % Energy-based approach for velocity calculation
    omega = zeros(size(theta));
    omega(1) = omega_avg;

    for i = 2:length(theta)
        dtheta = theta(i) - theta(i-1);
        dW = 0.5 * (T_net(i) + T_net(i-1)) * dtheta;
        omega_squared = omega(i-1)^2 + 2 * dW / I_flywheel;

        if omega_squared > 0
            omega(i) = sqrt(omega_squared);
        else
            omega(i) = 0.1 * omega_avg;
        end
    end

    % Adjust to maintain correct average speed
    omega_actual_avg = mean(omega);
    omega = omega * (omega_avg / omega_actual_avg);

    rpm = omega * 60 / (2*pi);

    omega_max = max(omega);
    omega_min = min(omega);
    omega_mean = mean(omega);
    Cs_actual = (omega_max - omega_min) / omega_mean;

    if any(omega <= 0)
        error('Non-positive angular velocity detected');
    end
end