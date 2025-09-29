clear; clc; close all;

addpath('..');

params = engine_parameters();
theta = linspace(0, 2*pi, params.simulationPointsPerCycle)';

[V_total, V_exp, V_comp, ~, ~] = calc_volumes(theta, params);
[P, m_total, P_mean] = schmidt_analysis(theta, V_total, V_exp, V_comp, params);

fprintf('\n=== Clean Code Results ===\n');
fprintf('Pressure Range: %.2f - %.2f kPa\n', min(P)/1000, max(P)/1000);
fprintf('Total Volume Range: %.2f - %.2f cm³\n', min(V_total)*1e6, max(V_total)*1e6);
fprintf('Compression Ratio: %.2f\n', max(V_total)/min(V_total));
fprintf('Hot Volume Range: %.2f - %.2f cm³\n', min(V_exp)*1e6, max(V_exp)*1e6);
fprintf('Cold Volume Range: %.2f - %.2f cm³\n', min(V_comp)*1e6, max(V_comp)*1e6);
fprintf('Regenerator Volume: %.2f cm³ (constant)\n', params.regeneratorVolume*1e6);

fprintf('\n=== PeytonTest Results (for comparison) ===\n');
fprintf('Pressure Range: 474.63 - 1178.83 kPa\n');
fprintf('Total Volume Range: 140.25 - 238.42 cm³\n');
fprintf('Compression Ratio: 1.70\n');
fprintf('Hot Volume Range: 32.26 - 110.80 cm³\n');
fprintf('Cold Volume Range: 35.61 - 163.18 cm³\n');
fprintf('Regenerator Volume: 20.00 cm³ (constant)\n');

figure('Name', 'P-V Diagram Comparison', 'Position', [100, 100, 600, 500]);
plot(V_total*1e6, P/1000, 'b-', 'LineWidth', 2);
xlabel('Total Volume (cm³)');
ylabel('Pressure (kPa)');
title('P-V Diagram - Updated Clean Code');
grid on;
axis tight;