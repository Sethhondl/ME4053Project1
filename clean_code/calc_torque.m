function [T_total, T_power, T_disp, T_mean] = calc_torque(P, theta, x_power, x_disp, params)
    % Calculate instantaneous torque from pressure forces

    r_p = params.powerCrankLength;
    l_p = params.powerRodLength;
    r_d = params.displacerCrankLength;
    l_d = params.displacerRodLength;
    phase = params.phaseShift;
    A = params.cylinderArea;
    P_atm = params.atmosphericPressure;

    T_power = zeros(size(theta));
    T_disp = zeros(size(theta));

    for i = 1:length(theta)
        % Power piston torque
        F_power = (P(i) - P_atm) * A;
        sin_beta_p = r_p * sin(theta(i)) / l_p;
        if abs(sin_beta_p) < 1
            cos_beta_p = sqrt(1 - sin_beta_p^2);
            T_power(i) = F_power * r_p * sin(theta(i)) / cos_beta_p;
        else
            T_power(i) = 0;
        end

        % Displacer torque
        theta_disp = theta(i) - phase;
        d_rod = params.displacerRodDiameter;
        A_rod = pi * (d_rod/2)^2;
        F_disp = (P(i) - P_atm) * A_rod;
        sin_beta_d = r_d * sin(theta_disp) / l_d;
        if abs(sin_beta_d) < 1
            cos_beta_d = sqrt(1 - sin_beta_d^2);
            T_disp(i) = F_disp * r_d * sin(theta_disp) / cos_beta_d;
        else
            T_disp(i) = 0;
        end
    end

    T_total = T_power + T_disp;
    T_mean = mean(T_total);

    if all(T_total == 0)
        error('Zero torque throughout cycle');
    end
end