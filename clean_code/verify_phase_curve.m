% Verify the phase optimization curve shape
params = engine_parameters();
[optimal_phase, energy_curve, power_curve, efficiency_curve] = optimize_phase(params);

% Check the power variation
fprintf('Power curve analysis:\n');
fprintf('  Number of points: %d\n', length(power_curve));
fprintf('  Min power: %.1f W\n', min(power_curve));
fprintf('  Max power: %.1f W\n', max(power_curve));
fprintf('  Power range: %.1f W\n', max(power_curve) - min(power_curve));
fprintf('  Relative variation: %.1f%%\n', 100*(max(power_curve)-min(power_curve))/mean(power_curve));

% Create the corrected phase range
coarse_phases = 60:5:120;
fine_phases = 100:0.1:110;
phase_range = unique([coarse_phases, fine_phases]);

% Plot to verify
figure('Name', 'Corrected Phase Optimization');
yyaxis left;
plot(phase_range, power_curve/1000, 'b-', 'LineWidth', 2);
ylabel('Power Output (kW)');
ylim([0.18 0.27]);
yyaxis right;
plot(phase_range, energy_curve, 'r--', 'LineWidth', 2);
ylabel('Work per Cycle (J)');
ylim([16 25]);
xlabel('Phase Angle (degrees)');
title('Phase Angle Optimization - Corrected');
grid on;
xlim([55, 125]);
xline(optimal_phase, 'k:', 'LineWidth', 1.5);
legend('Power Output', 'Work per Cycle', sprintf('Optimal: %.3f°', optimal_phase), 'Location', 'south');

fprintf('\nThe curve should now show:\n');
fprintf('  - Clear rise from 60° to ~105°\n');
fprintf('  - Peak around 105°\n');
fprintf('  - Decline from 105° to 120°\n');
fprintf('  - No artificial flat top\n');