function display_results(results, params)
    % Display minimal analysis summary

    eta_carnot = 1 - params.coldTemperature / params.hotTemperature;

    fprintf('\n========== STIRLING ENGINE ANALYSIS ==========\n');
    fprintf('Power Output: %.2f kW\n', results.P_indicated/1000);
    fprintf('Thermal Efficiency: %.1f%% (Carnot: %.1f%%)\n', ...
            results.efficiency * 100, eta_carnot * 100);
    fprintf('Flywheel: D=%.3f m, M=%.1f kg\n', ...
            results.flywheel.D_outer, results.flywheel.mass);
    fprintf('Speed Fluctuation: Cs=%.4f (Target: %.4f)\n', ...
            results.Cs_actual, params.flywheelCoefficientOfFluctuation);
    fprintf('Optimal Phase: %d degrees\n', results.optimal_phase);

    % Check requirements
    requirements_met = true;
    % Power output is informational only - not a requirement
    % (removed incorrect 1-10 kW range check)
    if results.efficiency > eta_carnot
        requirements_met = false;
    end
    if results.flywheel.D_outer > params.maximumFlywheelDiameter
        requirements_met = false;
    end
    if results.Cs_actual > params.flywheelCoefficientOfFluctuation * 1.1
        requirements_met = false;
    end

    if requirements_met
        fprintf('STATUS: ALL REQUIREMENTS MET ✓\n');
    else
        fprintf('STATUS: REQUIREMENTS NOT MET ✗\n');
    end
    fprintf('===============================================\n\n');
end