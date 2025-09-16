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
    
    % Ideal Stirling cycle - properly calculated with mass conservation
    % Process 1-2: Isothermal compression at T_cold
    % Process 2-3: Isochoric (constant volume) heating
    % Process 3-4: Isothermal expansion at T_hot
    % Process 4-1: Isochoric (constant volume) cooling

    % Use the same total mass as the actual engine for fair comparison
    m_total_ideal = results.m_total;  % Use actual engine's mass
    R = params.gasConstant;

    % Calculate ideal cycle pressures to match actual cycle's operating range
    % State 1: Beginning of compression (V_max, T_cold)
    % Use ideal gas law: P1 = m*R*T_cold/V_max
    % But we need to determine the mass distribution

    % For ideal Stirling cycle with perfect heat transfer:
    % At state 1: All gas at T_cold, volume V_max
    P1 = m_total_ideal * R * params.coldTemperature / V_max;

    % State 2: End of compression (V_min, T_cold)
    % Isothermal compression at T_cold: P1*V_max = P2*V_min
    P2 = P1 * (V_max / V_min);

    % State 3: After isochoric heating (V_min, T_hot)
    % Isochoric process: P2/T_cold = P3/T_hot
    P3 = P2 * (params.hotTemperature / params.coldTemperature);

    % State 4: After isothermal expansion (V_max, T_hot)
    % Isothermal expansion at T_hot: P3*V_min = P4*V_max
    P4 = P3 * (V_min / V_max);

    % Verify closure: P4/T_hot should equal P1/T_cold
    % P4 = P3*(V_min/V_max) = P2*(T_hot/T_cold)*(V_min/V_max)
    %    = P1*(V_max/V_min)*(T_hot/T_cold)*(V_min/V_max) = P1*(T_hot/T_cold)
    % So P4/T_hot = P1/T_cold ✓

    % Create arrays for plotting with more points for isothermal processes
    n_iso = 50;  % Points for isothermal curves

    % Process 1-2: Isothermal compression
    V_12 = linspace(V_max, V_min, n_iso);
    P_12 = P1 * V_max ./ V_12;

    % Process 2-3: Isochoric heating (vertical line)
    % Create more points for smoother visualization
    n_vert = 20;
    V_23 = ones(1, n_vert) * V_min;
    P_23 = linspace(P2, P3, n_vert);

    % Process 3-4: Isothermal expansion
    V_34 = linspace(V_min, V_max, n_iso);
    P_34 = P3 * V_min ./ V_34;

    % Process 4-1: Isochoric cooling (vertical line)
    % Create more points for smoother visualization
    V_41 = ones(1, n_vert) * V_max;
    P_41 = linspace(P4, P1, n_vert);

    % Combine for complete cycle
    % Add first point again at end to close the cycle
    V_ideal = [V_12, V_23, V_34, V_41, V_12(1)];
    P_ideal = [P_12, P_23, P_34, P_41, P_12(1)];
    
    % Plot actual engine cycle
    plot(results.V_total*1000, results.P/1e6, 'b-', 'LineWidth', 2);
    hold on;

    % Plot ideal Stirling cycle
    plot(V_ideal*1000, P_ideal/1e6, 'r--', 'LineWidth', 1.5);

    % Note: Direction arrows removed to avoid visual confusion
    % The cycle direction is clockwise (power-producing)

    % Add process labels with better positioning
    text(mean([V_min, V_max])*1000, (P1/1e6 + P2/1e6)/2, ...
         '1→2: Compression (T_{cold})', ...
         'HorizontalAlignment', 'center', 'FontSize', 9);
    text(V_min*1000*0.97, mean([P2, P3])/1e6, '2→3', ...
         'HorizontalAlignment', 'right', 'FontSize', 9);
    text(mean([V_min, V_max])*1000, (P3/1e6 + P4/1e6)/2, ...
         '3→4: Expansion (T_{hot})', ...
         'HorizontalAlignment', 'center', 'FontSize', 9);
    text(V_max*1000*1.03, mean([P4, P1])/1e6, '4→1', ...
         'HorizontalAlignment', 'left', 'FontSize', 9);
    
    % Format plot
    xlabel('Volume (L)', 'FontSize', 12);
    ylabel('Pressure (MPa)', 'FontSize', 12);
    title('P-V Diagram: Stirling Cycle vs Engine Cycle', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Engine Cycle (Schmidt)', 'Ideal Stirling Cycle', 'Location', 'best');
    grid on;

    % Set axis limits to show only positive values
    xlim([0, max(results.V_total)*1000*1.1]);  % 0 to max volume plus 10% margin
    ylim([0, max([max(results.P), max(P_ideal)])/1e6*1.1]);  % 0 to max pressure plus 10% margin
    
    % Add work area annotation
    work_actual = abs(results.W_indicated);
    % Calculate ideal Stirling cycle work
    % W = m*R*(T_hot - T_cold)*ln(V_max/V_min)
    if isfield(results, 'm_total')
        work_ideal = results.m_total * params.gasConstant * (params.hotTemperature - params.coldTemperature) * log(V_max/V_min);
    else
        work_ideal = work_actual * 1.5;  % Estimate
    end

    % Position text in empty area of plot
    text_x = V_min*1000 + (V_max - V_min)*1000*0.7;
    text_y = P_min/1e6 + (P_max - P_min)/1e6*0.8;
    text(text_x, text_y, ...
         sprintf('Work per Cycle:\nActual = %.1f J\nIdeal = %.1f J', work_actual, work_ideal), ...
         'FontSize', 9, 'BackgroundColor', 'white', 'EdgeColor', 'black');

    % Add cycle direction indicator in a visible location
    text(V_max*1000*0.9, P_max/1e6*0.9, 'Clockwise = Power', ...
         'FontSize', 9, 'FontWeight', 'bold', 'Color', 'blue');
    
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
    current_phase = params.phaseShift * 180/pi;
    plot(current_phase, interp1(results.optimization.energy_curve(:,1), ...
         results.optimization.energy_curve(:,2), current_phase), ...
         'gs', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'LineWidth', 2);
    
    xlabel('Phase Angle (degrees)', 'FontSize', 11);
    ylabel('Energy per Cycle (J)', 'FontSize', 11);
    title('Energy and Power vs Phase Angle', 'FontSize', 13, 'FontWeight', 'bold');
    legend('Energy', sprintf('Optimal (%.0f°)', results.optimal_phase), ...
           sprintf('Current (%.0f°)', current_phase), 'Location', 'best');
    grid on;
    xlim([45, 135]);  % Expanded range
    
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
    xlim([45, 135]);  % Expanded range
    
    % Save plot
    saveas(gcf, 'results/phase_optimization.png');
    
    fprintf('All plots have been generated and saved to the results/ directory\n');
    fprintf('Files created:\n');
    fprintf('  - pv_diagram.png\n');
    fprintf('  - torque_profile.png\n');
    fprintf('  - velocity_variation.png\n');
    fprintf('  - phase_optimization.png\n\n');
end