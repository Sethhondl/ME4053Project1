params = engine_parameters();
theta = linspace(0, 2*pi, 361);
[x_power, x_disp, V_total, V_exp, V_comp] = calc_volumes(theta, params);
[P, m_total, P_mean] = schmidt_analysis(theta, V_total, V_exp, V_comp, params);
[T_total, T_power, T_disp, T_mean] = calc_torque(P, theta, x_power, x_disp, params);

fprintf('Mean torque: %.2f N·m\n', T_mean);
fprintf('Mean positive torque: %.2f N·m\n', mean(T_total(T_total > 0)));
fprintf('Mean negative torque: %.2f N·m\n', mean(T_total(T_total < 0)));
fprintf('Positive torque samples: %d\n', sum(T_total > 0));
fprintf('Negative torque samples: %d\n', sum(T_total < 0));

% Check power calculation
W_pdv = trapz(theta, T_total) * params.averageRPM * 2 * pi / 60;
fprintf('\nPower from torque: %.2f W\n', W_pdv);