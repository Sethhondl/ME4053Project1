% Test phase optimization plot
params = engine_parameters();
[optimal_phase, energy_curve, power_curve, efficiency_curve] = optimize_phase(params);

fprintf('Graph verification:\n');
fprintf('  Optimal phase: %.3f degrees\n', optimal_phase);
fprintf('  Power curve has %d points\n', length(power_curve));
fprintf('  X-axis should show: 60 to 120 degrees\n');
fprintf('  Vertical line should be at: %.3f degrees (NOT %.0f)\n', optimal_phase, optimal_phase * 180/pi);
fprintf('  Legend should show: Optimal: %.3f°\n', optimal_phase);

% Create a test plot to verify
figure('Name', 'Test Phase Plot');
phase_range = linspace(60, 120, length(power_curve));
yyaxis left;
plot(phase_range, power_curve/1000, 'b-', 'LineWidth', 2);
ylabel('Power Output (kW)');
yyaxis right;
plot(phase_range, energy_curve, 'r--', 'LineWidth', 2);
ylabel('Work per Cycle (J)');
xlabel('Phase Angle (degrees)');
title('Phase Angle Optimization - TEST');
grid on;
xlim([55, 125]);  % Fixed x-axis limits
xline(optimal_phase, 'k:', 'LineWidth', 1.5);  % Correct position
legend('Power Output', 'Work per Cycle', sprintf('Optimal: %.3f°', optimal_phase), 'Location', 'south');

fprintf('\nIf the plot looks correct with x-axis from 55-125 degrees\n');
fprintf('and vertical line at ~105 degrees, the fix is working!\n');