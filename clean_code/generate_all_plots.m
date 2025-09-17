function generate_all_plots(results, params)
    % Generate all required plots

    if ~exist('results', 'dir')
        mkdir('results');
    end

    % P-V Diagram
    figure('Position', [100, 100, 800, 600]);
    V_min = min(results.V_total);
    V_max = max(results.V_total);

    % Ideal Stirling cycle
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

    V_ideal = [V_12, V_23, V_34, V_41, V_12(1)];
    P_ideal = [P_12, P_23, P_34, P_41, P_12(1)];

    plot(results.V_total*1000, results.P/1e6, 'b-', 'LineWidth', 2);
    hold on;
    plot(V_ideal*1000, P_ideal/1e6, 'r--', 'LineWidth', 1.5);
    xlabel('Volume (L)', 'FontSize', 12);
    ylabel('Pressure (MPa)', 'FontSize', 12);
    title('P-V Diagram', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Engine Cycle', 'Ideal Stirling', 'Location', 'best');
    grid on;
    xlim([0, max(results.V_total)*1000*1.1]);
    ylim([0, max([max(results.P), max(P_ideal)])/1e6*1.1]);
    saveas(gcf, 'results/pv_diagram.png');

    % Torque vs Crank Angle
    figure('Position', [200, 150, 800, 600]);
    theta_deg = results.theta * 180/pi;
    plot(theta_deg, results.T_total, 'b-', 'LineWidth', 2);
    hold on;
    plot(theta_deg, results.T_power, 'g--', 'LineWidth', 1);
    plot(theta_deg, results.T_disp, 'r--', 'LineWidth', 1);
    plot([0, 360], [results.T_mean, results.T_mean], 'k:', 'LineWidth', 1.5);
    xlabel('Crank Angle (degrees)', 'FontSize', 12);
    ylabel('Torque (N·m)', 'FontSize', 12);
    title('Torque Profile', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Total', 'Power Piston', 'Displacer', 'Mean', 'Location', 'best');
    grid on;
    xlim([0, 360]);
    saveas(gcf, 'results/torque_profile.png');

    % Velocity Variation
    figure('Position', [300, 200, 800, 600]);
    plot(theta_deg, results.rpm, 'b-', 'LineWidth', 2);
    hold on;
    plot([0, 360], [results.rpm_mean, results.rpm_mean], 'r--', 'LineWidth', 1.5);
    xlabel('Crank Angle (degrees)', 'FontSize', 12);
    ylabel('Speed (RPM)', 'FontSize', 12);
    title('Speed Variation', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Instantaneous', 'Mean', 'Location', 'best');
    grid on;
    xlim([0, 360]);
    saveas(gcf, 'results/velocity_variation.png');

    % Phase Optimization
    figure('Position', [400, 250, 800, 600]);
    subplot(2,1,1);
    plot(results.optimization.energy_curve(:,1), ...
        results.optimization.energy_curve(:,2), 'b-o', 'LineWidth', 2);
    xlabel('Phase Angle (degrees)', 'FontSize', 11);
    ylabel('Energy per Cycle (J)', 'FontSize', 11);
    title('Phase Optimization', 'FontSize', 13, 'FontWeight', 'bold');
    grid on;
    xlim([45, 135]);

    subplot(2,1,2);
    plot(results.optimization.power_curve(:,1), ...
        results.optimization.power_curve(:,2)/1000, 'r-o', 'LineWidth', 2);
    xlabel('Phase Angle (degrees)', 'FontSize', 11);
    ylabel('Power Output (kW)', 'FontSize', 11);
    grid on;
    xlim([45, 135]);
    saveas(gcf, 'results/phase_optimization.png');
end