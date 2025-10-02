function [D_outer, I_required, mass, energy_fluctuation] = size_flywheel(T_total, theta, params)
    omega_avg = params.averageRPM * 2 * pi / 60;
    Cs = params.flywheelCoefficientOfFluctuation;
    w = params.flywheelWidth;
    t = params.flywheelRimThickness;
    rho = params.flywheelMaterialDensity;

    T_mean = mean(T_total);
    T_deviation = T_total - T_mean;

    energy_variation = zeros(size(theta));
    for i = 2:length(theta)
        dtheta = theta(i) - theta(i-1);
        energy_variation(i) = energy_variation(i-1) + 0.5 * (T_deviation(i) + T_deviation(i-1)) * dtheta;
    end

    E_max = max(energy_variation);
    E_min = min(energy_variation);
    energy_fluctuation = E_max - E_min;

    I_initial = energy_fluctuation / (Cs * omega_avg^2);
    I_required = I_initial;

    for iter_cs = 1:10
        [~, ~, ~, Cs_actual] = simulate_dynamics(T_total, theta, I_required, params);

        error_cs = abs(Cs_actual - Cs) / Cs;
        if error_cs < 0.0001
            break;
        end

        correction_factor = Cs_actual / Cs;
        I_required = I_required * correction_factor;
    end

    r_outer_guess = (I_required / (2 * pi * rho * w * t))^(1/3);

    for iter = 1:10
        r_inner = r_outer_guess - t;
        if r_inner <= 0
            error('Flywheel thickness too large for required inertia');
        end

        V = pi * w * (r_outer_guess^2 - r_inner^2);
        m = rho * V;
        I_calc = 0.5 * m * (r_outer_guess^2 + r_inner^2);

        error_ratio = I_required / I_calc;
        r_outer_guess = r_outer_guess * error_ratio^(1/3);

        if abs(I_calc - I_required) / I_required < 0.001
            break;
        end

    end

    r_outer = r_outer_guess;
    r_inner = r_outer - t;
    V = pi * w * (r_outer^2 - r_inner^2);
    mass = rho * V;
    D_outer = 2 * r_outer;

    if D_outer > params.maximumFlywheelDiameter
        warning('Flywheel diameter exceeds practical limit of %.1f m', params.maximumFlywheelDiameter);
    end
end
end