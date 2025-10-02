params = engine_parameters();
theta = linspace(0, 2*pi, 361);
[x_power, x_disp, V_total, V_exp, V_comp] = calc_volumes(theta, params);
[P, m_total, P_mean] = schmidt_analysis(theta, V_total, V_exp, V_comp, params);
[T_total, T_power, T_disp, T_mean] = calc_torque(P, theta, x_power, x_disp, params);

% Power from mean torque
omega = params.averageRPM * 2 * pi / 60;  % rad/s
P_from_mean_torque = abs(T_mean) * omega;
fprintf('Power from mean torque: %.2f W\n', P_from_mean_torque);

% Power from torque integration over cycle
W_cycle = trapz(theta, T_total);  % Work per cycle (J)
f = params.averageRPM / 60;  % Hz
P_from_integration = abs(W_cycle) * f;
fprintf('Power from torque integration: %.2f W\n', P_from_integration);

% Calculate actual power from P-V
[W_indicated, P_indicated, W_mep, P_mep, MEP, efficiency] = calc_power(P, V_total, theta, params);
fprintf('Power from P-V (Method 1): %.2f W\n', P_indicated);
fprintf('Power from P-V (Method 2): %.2f W\n', P_mep);