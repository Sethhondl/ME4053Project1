% Debug ideal cycle
params = engine_parameters();
V_min = 0.0002;
V_max = 0.003254;
P_min_actual = 0.083e6;
P_max_actual = 1.5e6;

% Calculate ideal cycle pressures
P1 = P_min_actual * 1.2;  % Start slightly above minimum

% Isothermal compression
P2 = P1 * (V_max / V_min);
fprintf('P1 = %.3f MPa, P2 = %.3f MPa\n', P1/1e6, P2/1e6);

% Isochoric heating
P3 = P2 * (params.T_hot / params.T_cold);
fprintf('P3 = %.3f MPa\n', P3/1e6);

% Isothermal expansion
P4 = P3 * (V_min / V_max);
fprintf('P4 = %.3f MPa\n', P4/1e6);

% Check closure
P1_check = P4 * (params.T_cold / params.T_hot);
fprintf('P1 check = %.3f MPa (should equal P1)\n', P1_check/1e6);

% All pressures should be positive
fprintf('\nAll positive? P1>0: %d, P2>0: %d, P3>0: %d, P4>0: %d\n', ...
    P1>0, P2>0, P3>0, P4>0);
