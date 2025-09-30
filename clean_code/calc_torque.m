function [T_total, T_power, T_disp, T_mean] = calc_torque(P, theta, x_power, x_disp, params)
    r_p = params.powerCrankLength;
    l_p = params.powerRodLength;
    r_d = params.displacerCrankLength;
    l_d = params.displacerRodLength;
    phase = params.phaseShift;
    A = pi * (params.cylinderBore/2)^2;
    P_atm = params.atmosphericPressure;

    T_power = zeros(size(theta));
    T_disp = zeros(size(theta));

    for i = 1:length(theta)
        F_power = (P(i) - P_atm) * A;

        sin_beta_p = r_p * sin(theta(i)) / l_p;
        if abs(sin_beta_p) < 1
            cos_beta_p = sqrt(1 - sin_beta_p^2);
            T_power(i) = -F_power * r_p * sin(theta(i)) / cos_beta_p;
        else
            T_power(i) = 0;
        end

        T_disp(i) = 0;
    end

    T_total = T_power + T_disp;
    T_mean = mean(T_total);
end