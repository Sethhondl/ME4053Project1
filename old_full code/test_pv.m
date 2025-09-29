% Test ideal Stirling cycle calculation
V_min = 0.0002;
V_max = 0.003254;
P_mean = 0.79e6;
T_hot = 600;
T_cold = 350;

P1 = P_mean * sqrt(V_min/V_max);
P2 = P1 * (V_max/V_min);
P3 = P2 * (T_hot/T_cold);
P4 = P3 * (V_min/V_max);

fprintf('Ideal Stirling Cycle Pressures:\n');
fprintf('P1 = %.3f MPa\n', P1/1e6);
fprintf('P2 = %.3f MPa\n', P2/1e6);
fprintf('P3 = %.3f MPa\n', P3/1e6);
fprintf('P4 = %.3f MPa\n', P4/1e6);
fprintf('Compression ratio = %.2f\n', V_max/V_min);
