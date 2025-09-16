function generate_all_plots(results, params)
    % GENERATE_ALL_PLOTS Create all required plots for the Stirling engine analysis
    %
    % Inputs:
    %   results - structure containing all analysis results
    %   params - engine parameters structure
    
    % Create results directory if it doesn't exist
    if ~exist('results', 'dir')
        mkdir('results');
    end
    
    %% Plot 1: P-V Diagram (Stirling Cycle vs Engine Cycle)
    figure('Position', [100, 100, 800, 600]);
    
    % Calculate ideal Stirling cycle for comparison
    V_min = min(results.V_total);
    V_max = max(results.V_total);
    P_min = min(results.P);
    P_max = max(results.P);
    
    % Ideal Stirling cycle points (rectangular on P-V diagram)
    % Process 1-2: Isothermal compression at T_cold
    % Process 2-3: Constant volume heating
    % Process 3-4: Isothermal expansion at T_hot
    % Process 4-1: Constant volume cooling
    
    % Define ideal cycle corners
    V_ideal = [V_min, V_min, V_max, V_max, V_min];
    
    % For ideal isothermal processes: P*V = constant
    % At hot temperature: P3*V_min = P4*V_max
    % At cold temperature: P1*V_max = P2*V_min
    
    % Approximate ideal pressures based on temperature ratio
    T_ratio = params.T_hot / params.T_cold;
    P2 = P_max;  % Maximum pressure at minimum volume, cold temperature
    P3 = P2 * T_ratio;  % Pressure after constant volume heating
    P4 = P3 * (V_min/V_max);  % Pressure after isothermal expansion
    P1 = P4 / T_ratio;  % Pressure after constant volume cooling
    
    P_ideal = [P1, P2, P3, P4, P1];
    
    % Plot actual engine cycle
    plot(results.V_total*1000, results.P/1e6, 'b-', 'LineWidth', 2);
    hold on;
    
    % Plot ideal Stirling cycle
    plot(V_ideal*1000, P_ideal/1e6, 'r--', 'LineWidth', 1.5);
    
    % Add process labels
    text(mean([V_min, V_max])*1000, P1/1e6*0.9, '1: Isothermal Compression', ...
         'HorizontalAlignment', 'center', 'FontSize', 9);
    text(V_min*1000*0.95, mean([P2, P3])/1e6, '2: Heating', ...
         'Rotation', 90, 'HorizontalAlignment', 'center', 'FontSize', 9);
    text(mean([V_min, V_max])*1000, P3/1e6*1.1, '3: Isothermal Expansion', ...
         'HorizontalAlignment', 'center', 'FontSize', 9);
    text(V_max*1000*1.05, mean([P4, P1])/1e6, '4: Cooling', ...
         'Rotation', 90, 'HorizontalAlignment', 'center', 'FontSize', 9);
    
    % Format plot
    xlabel('Volume (L)', 'FontSize', 12);
    ylabel('Pressure (MPa)', 'FontSize', 12);
    title('P-V Diagram: Stirling Cycle vs Engine Cycle', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Engine Cycle (Schmidt)', 'Ideal Stirling Cycle', 'Location', 'best');
    grid on;
    
    % Add work area annotation
    work_area = abs(trapz(results.V_total, results.P));
    text(V_min*1000*1.1, P_min/1e6*1.2, ...
         sprintf('Work = %.1f J', work_area), 'FontSize', 10);
    
    % Save plot
    saveas(gcf, 'results/pv_diagram.png');
    
    %% Plot 2: Torque vs Crank Angle
    figure('Position', [200, 150, 800, 600]);
    
    theta_deg = results.theta * 180/pi;
    plot(theta_deg, results.T_total, 'b-', 'LineWidth', 2);
    hold on;
    plot(theta_deg, results.T_power, 'g--', 'LineWidth', 1);
    plot(theta_deg, results.T_disp, 'r--', 'LineWidth', 1);
    plot([0, 360], [results.T_mean, results.T_mean], 'k:', 'LineWidth', 1.5);
    
    % Mark zero crossings
    zero_crossings = find(diff(sign(results.T_total)) ~= 0);
    if ~isempty(zero_crossings)
        plot(theta_deg(zero_crossings), results.T_total(zero_crossings), ...
             'ro', 'MarkerSize', 8, 'LineWidth', 2);
    end
    
    xlabel('Crank Angle (degrees)', 'FontSize', 12);
    ylabel('Torque (N·m)', 'FontSize', 12);
    title('Torque vs Crank Angle', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Total Torque', 'Power Piston', 'Displacer', 'Mean Torque', ...
           'Location', 'best');
    grid on;
    xlim([0, 360]);
    
    % Add annotations
    text(300, results.T_mean*1.2, sprintf('Mean = %.1f N·m', results.T_mean), ...
         'FontSize', 10);
    
    % Save plot
    saveas(gcf, 'results/torque_profile.png');
    
    %% Plot 3: Rotational Velocity vs Crank Angle
    figure('Position', [300, 200, 800, 600]);
    
    plot(theta_deg, results.rpm, 'b-', 'LineWidth', 2);
    hold on;
    plot([0, 360], [results.rpm_mean, results.rpm_mean], 'r--', 'LineWidth', 1.5);
    
    % Mark max and min speeds
    [rpm_max, max_idx] = max(results.rpm);
    [rpm_min, min_idx] = min(results.rpm);
    plot(theta_deg(max_idx), rpm_max, 'g^', 'MarkerSize', 10, ...
         'MarkerFaceColor', 'g', 'LineWidth', 2);
    plot(theta_deg(min_idx), rpm_min, 'rv', 'MarkerSize', 10, ...
         'MarkerFaceColor', 'r', 'LineWidth', 2);
    
    xlabel('Crank Angle (degrees)', 'FontSize', 12);
    ylabel('Rotational Speed (RPM)', 'FontSize', 12);
    title('Rotational Velocity vs Crank Angle', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Instantaneous Speed', 'Mean Speed', 'Maximum', 'Minimum', ...
           'Location', 'best');
    grid on;
    xlim([0, 360]);
    
    % Add coefficient of fluctuation annotation
    text(180, rpm_min*0.95, ...
         sprintf('Cs = %.4f', results.Cs_actual), ...
         'FontSize', 11, 'HorizontalAlignment', 'center', ...
         'BackgroundColor', 'white', 'EdgeColor', 'black');
    
    % Add speed variation info
    speed_var = rpm_max - rpm_min;
    text(270, rpm_max*0.98, ...
         sprintf('ΔN = %.1f RPM', speed_var), ...
         'FontSize', 10);
    
    % Save plot
    saveas(gcf, 'results/velocity_variation.png');
    
    %% Plot 4: Energy per Cycle vs Phase Angle
    figure('Position', [400, 250, 800, 600]);
    
    % Create subplot for energy and power
    subplot(2,1,1);
    plot(results.optimization.energy_curve(:,1), ...
         results.optimization.energy_curve(:,2), 'b-o', 'LineWidth', 2);
    hold on;
    
    % Mark optimal point
    [~, opt_idx] = max(results.optimization.energy_curve(:,2));
    plot(results.optimization.energy_curve(opt_idx,1), ...
         results.optimization.energy_curve(opt_idx,2), ...
         'r*', 'MarkerSize', 15, 'LineWidth', 2);
    
    % Mark current setting
    current_phase = params.phase_shift * 180/pi;
    plot(current_phase, interp1(results.optimization.energy_curve(:,1), ...
         results.optimization.energy_curve(:,2), current_phase), ...
         'gs', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'LineWidth', 2);
    
    xlabel('Phase Angle (degrees)', 'FontSize', 11);
    ylabel('Energy per Cycle (J)', 'FontSize', 11);
    title('Energy and Power vs Phase Angle', 'FontSize', 13, 'FontWeight', 'bold');
    legend('Energy', sprintf('Optimal (%.0f°)', results.optimal_phase), ...
           sprintf('Current (%.0f°)', current_phase), 'Location', 'best');
    grid on;
    xlim([60, 120]);
    
    % Power subplot
    subplot(2,1,2);
    plot(results.optimization.power_curve(:,1), ...
         results.optimization.power_curve(:,2)/1000, 'r-o', 'LineWidth', 2);
    hold on;
    
    % Mark optimal point
    [~, opt_idx] = max(results.optimization.power_curve(:,2));
    plot(results.optimization.power_curve(opt_idx,1), ...
         results.optimization.power_curve(opt_idx,2)/1000, ...
         'b*', 'MarkerSize', 15, 'LineWidth', 2);
    
    xlabel('Phase Angle (degrees)', 'FontSize', 11);
    ylabel('Power Output (kW)', 'FontSize', 11);
    legend('Power', sprintf('Maximum (%.0f°)', ...
           results.optimization.power_curve(opt_idx,1)), 'Location', 'best');
    grid on;
    xlim([60, 120]);
    
    % Save plot
    saveas(gcf, 'results/phase_optimization.png');
    
    fprintf('All plots have been generated and saved to the results/ directory\n');
    fprintf('Files created:\n');
    fprintf('  - pv_diagram.png\n');
    fprintf('  - torque_profile.png\n');
    fprintf('  - velocity_variation.png\n');
    fprintf('  - phase_optimization.png\n\n');
end