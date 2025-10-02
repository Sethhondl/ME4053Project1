params = engine_parameters();
theta = linspace(0, 2*pi, 361);
[x_power, x_disp, V_total, V_exp, V_comp] = calc_volumes(theta, params);
[P, m_total, P_mean] = schmidt_analysis(theta, V_total, V_exp, V_comp, params);

% Calculate work from P-V diagram
W_pv = -trapz(V_total, P);  % Negative because we integrate P dV
fprintf('Work from P-V diagram: %.2f J\n', W_pv);

% Check if cycle is clockwise (power producing)
mid = floor(length(theta)/2);
area_test = trapz(V_total(1:mid), P(1:mid)) - trapz(V_total(mid:end), P(mid:end));
if area_test < 0
    fprintf('Cycle is clockwise (power producing)\n');
else
    fprintf('Cycle is counter-clockwise (power consuming)\n');
end

% Power from P-V work
operatingFrequency = params.averageRPM / 60;
P_from_pv = W_pv * operatingFrequency;
fprintf('Power from P-V work: %.2f W\n', P_from_pv);