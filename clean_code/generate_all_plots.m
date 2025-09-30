function generate_all_plots(results, params)
    if ~exist('results', 'dir')
        mkdir('results');
    end
    
    figure('Position', [100, 100, 800, 600]);
    V_min = min(results.V_total);
    V_max = max(results.V_total);
    m_total_ideal = results.m_total;
    R = params.gasConstant;

    P1 = m_total_ideal * R * params.coldTemperature / V_max;
    P2 = P1 * (V_max / V_min);
    P3 = P2 * (params.hotTemperature / params.coldTemperature);
    P4 = P3 * (V_min / V_max);

    n_iso = 50;
    V_12 = linspace(V_max, V_min, n_iso);
    P_12 = P1 * V_max ./ V_12;
    V_23 = ones(1, 20) * V_min;
    P_23 = linspace(P2, P3, 20);
    V_34 = linspace(V_min, V_max, n_iso);
    P_34 = P3 * V_min ./ V_34;
    V_41 = ones(1, 20) * V_max;
    P_41 = linspace(P4, P1, 20);

    V_ideal = [V_12, V_23, V_34, V_41];
    P_ideal = [P_12, P_23, P_34, P_41];

    plot(V_ideal * 1e6, P_ideal / 1e6, 'k--', 'LineWidth', 2);
    hold on;
    plot(results.V_total * 1e6, results.P / 1e6, 'b-', 'LineWidth', 2);
    xlabel('Volume (cm³)');
    ylabel('Pressure (MPa)');
    title('P-V Diagram');
    legend('Ideal Stirling Cycle', 'Actual Engine Cycle', 'Location', 'northeast');
    grid on;
    box on;
    saveas(gcf, 'results/pv_diagram.png');

    figure('Position', [100, 300, 800, 600]);
    plot(results.theta * 180/pi, results.T_total, 'b-', 'LineWidth', 2);
    hold on;
    plot(results.theta * 180/pi, ones(size(results.theta)) * results.T_mean, 'r--', 'LineWidth', 1.5);
    xlabel('Crank Angle (degrees)');
    ylabel('Torque (N·m)');
    title('Torque Profile');
    legend('Instantaneous Torque', 'Mean Torque', 'Location', 'northeast');
    grid on;
    xlim([0 360]);
    saveas(gcf, 'results/torque_profile.png');

    figure('Position', [100, 500, 800, 600]);
    % Create phase range matching the optimization data
    % Use combined coarse (60:5:120) and fine (100:0.1:110) resolution
    coarse_phases = 60:5:120;
    fine_phases = 100:0.1:110;
    phase_range = unique([coarse_phases, fine_phases]);

    % The optimization returns data at these phase points
    % Ensure we have the right number of points
    if length(phase_range) ~= length(results.optimization.power_curve)
        % If mismatch, use simple linear spacing as fallback
        phase_range = linspace(60, 120, length(results.optimization.power_curve));
    end

    yyaxis left;
    plot(phase_range, results.optimization.power_curve/1000, 'b-', 'LineWidth', 2);
    ylabel('Power Output (kW)');
    yyaxis right;
    plot(phase_range, results.optimization.energy_curve, 'r--', 'LineWidth', 2);
    ylabel('Work per Cycle (J)');
    xlabel('Phase Angle (degrees)');
    title('Phase Angle Optimization');
    grid on;
    xlim([55, 125]);
    xline(results.optimal_phase, 'k:', 'LineWidth', 1.5);
    legend('Power Output', 'Work per Cycle', sprintf('Optimal: %.3f°', results.optimal_phase), 'Location', 'south');
    saveas(gcf, 'results/phase_optimization.png');

    figure('Position', [100, 700, 800, 600]);
    plot(results.theta * 180/pi, results.rpm, 'b-', 'LineWidth', 2);
    hold on;
    plot(results.theta * 180/pi, ones(size(results.theta)) * results.rpm_mean, 'r--', 'LineWidth', 1.5);
    xlabel('Crank Angle (degrees)');
    ylabel('Speed (RPM)');
    title('Angular Velocity Variation');
    legend('Instantaneous Speed', 'Mean Speed', 'Location', 'northeast');
    grid on;
    xlim([0 360]);
    saveas(gcf, 'results/velocity_variation.png');
end
