% Test P-V diagram shapes with different dead volume ratios
clear all; close all;

% Load base parameters
params = engine_parameters();

% Save original dead volumes
original_dead_hot = params.deadVolumeHot;
original_dead_cold = params.deadVolumeCold;
original_regen = params.regeneratorVolume;

% Test cases with different dead volume ratios
figure('Position', [100, 100, 1200, 800]);

% Case 1: Current (CR = 1.7, large dead volume)
subplot(2,2,1);
theta = linspace(0, 2*pi, 360)';
[V_total1, V_exp1, V_comp1, ~, ~] = calc_volumes(theta, params);
[P1, ~, ~] = schmidt_analysis(theta, V_total1, V_exp1, V_comp1, params);
plot(V_total1*1000, P1/1e6, 'b-', 'LineWidth', 2);
xlabel('Volume (mL)'); ylabel('Pressure (MPa)');
title(sprintf('Current: CR=%.1f (Airfoil shape)', max(V_total1)/min(V_total1)));
grid on;
axis equal;

% Case 2: Reduced dead volume (CR ~ 2.5)
subplot(2,2,2);
params.deadVolumeHot = original_dead_hot * 0.3;
params.deadVolumeCold = original_dead_cold * 0.3;
params.regeneratorVolume = original_regen * 0.5;
params.totalDeadVolume = params.deadVolumeHot + params.deadVolumeCold + params.regeneratorVolume;
[V_total2, V_exp2, V_comp2, ~, ~] = calc_volumes(theta, params);
[P2, ~, ~] = schmidt_analysis(theta, V_total2, V_exp2, V_comp2, params);
plot(V_total2*1000, P2/1e6, 'r-', 'LineWidth', 2);
xlabel('Volume (mL)'); ylabel('Pressure (MPa)');
title(sprintf('Reduced dead: CR=%.1f', max(V_total2)/min(V_total2)));
grid on;
axis equal;

% Case 3: Minimal dead volume (CR ~ 4)
subplot(2,2,3);
params.deadVolumeHot = original_dead_hot * 0.1;
params.deadVolumeCold = original_dead_cold * 0.1;
params.regeneratorVolume = original_regen * 0.2;
params.totalDeadVolume = params.deadVolumeHot + params.deadVolumeCold + params.regeneratorVolume;
[V_total3, V_exp3, V_comp3, ~, ~] = calc_volumes(theta, params);
[P3, ~, ~] = schmidt_analysis(theta, V_total3, V_exp3, V_comp3, params);
plot(V_total3*1000, P3/1e6, 'g-', 'LineWidth', 2);
xlabel('Volume (mL)'); ylabel('Pressure (MPa)');
title(sprintf('Minimal dead: CR=%.1f', max(V_total3)/min(V_total3)));
grid on;
axis equal;

% Case 4: All three overlaid (normalized)
subplot(2,2,4);
% Normalize to show shapes clearly
V1_norm = (V_total1 - min(V_total1))/(max(V_total1) - min(V_total1));
P1_norm = (P1 - min(P1))/(max(P1) - min(P1));
V2_norm = (V_total2 - min(V_total2))/(max(V_total2) - min(V_total2));
P2_norm = (P2 - min(P2))/(max(P2) - min(P2));
V3_norm = (V_total3 - min(V_total3))/(max(V_total3) - min(V_total3));
P3_norm = (P3 - min(P3))/(max(P3) - min(P3));

plot(V1_norm, P1_norm, 'b-', 'LineWidth', 2); hold on;
plot(V2_norm, P2_norm, 'r-', 'LineWidth', 2);
plot(V3_norm, P3_norm, 'g-', 'LineWidth', 2);
xlabel('Normalized Volume'); ylabel('Normalized Pressure');
title('Shape Comparison (Normalized)');
legend('CR=1.7 (Current)', 'CR~2.5', 'CR~4', 'Location', 'best');
grid on;
axis equal;

sgtitle('P-V Diagram Shape vs Dead Volume Ratio');

% Print analysis
fprintf('\n=== P-V DIAGRAM SHAPE ANALYSIS ===\n');
fprintf('Case 1 (Current): CR = %.2f, Dead/Swept = %.2f\n', max(V_total1)/min(V_total1), params.totalDeadVolume/params.powerSweptVolume);
fprintf('  Shape: Narrow/Airfoil - Low work area\n');

% Reset for case 2 calculation
params.deadVolumeHot = original_dead_hot * 0.3;
params.deadVolumeCold = original_dead_cold * 0.3;
params.regeneratorVolume = original_regen * 0.5;
params.totalDeadVolume = params.deadVolumeHot + params.deadVolumeCold + params.regeneratorVolume;
fprintf('Case 2: CR = %.2f, Dead/Swept = %.2f\n', max(V_total2)/min(V_total2), params.totalDeadVolume/params.powerSweptVolume);
fprintf('  Shape: More rounded - Better work area\n');

% Reset for case 3 calculation
params.deadVolumeHot = original_dead_hot * 0.1;
params.deadVolumeCold = original_dead_cold * 0.1;
params.regeneratorVolume = original_regen * 0.2;
params.totalDeadVolume = params.deadVolumeHot + params.deadVolumeCold + params.regeneratorVolume;
fprintf('Case 3: CR = %.2f, Dead/Swept = %.2f\n', max(V_total3)/min(V_total3), params.totalDeadVolume/params.powerSweptVolume);
fprintf('  Shape: Kidney bean - Maximum work area\n');

% Calculate work for each case
W1 = -trapz(V_total1, P1);
W2 = -trapz(V_total2, P2);
W3 = -trapz(V_total3, P3);
fprintf('\nWork per cycle:\n');
fprintf('  Case 1: %.1f J\n', abs(W1));
fprintf('  Case 2: %.1f J\n', abs(W2));
fprintf('  Case 3: %.1f J\n', abs(W3));
fprintf('  Improvement: %.0f%% (Case 1 to 3)\n', (abs(W3)/abs(W1) - 1)*100);