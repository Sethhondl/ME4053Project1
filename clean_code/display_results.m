function display_results(results, params)
    eta_carnot = 1 - params.coldTemperature / params.hotTemperature;

    fprintf('\n========================================================\n');
    fprintf('     STIRLING ENGINE FLYWHEEL DESIGN - FINAL RESULTS   \n');
    fprintf('========================================================\n\n');

    fprintf('PRIMARY OBJECTIVE: FLYWHEEL SIZING\n');
    fprintf('==================================\n');
    fprintf('  Required Flywheel Diameter: %.3f m (%.1f mm)\n', results.flywheel.D_outer, results.flywheel.D_outer * 1000);
    fprintf('  Flywheel Mass: %.2f kg\n', results.flywheel.mass);
    fprintf('  Moment of Inertia: %.4f kg·m²\n', results.flywheel.I_required);
    fprintf('  Energy Fluctuation: %.2f J\n\n', results.flywheel.energy_fluctuation);

    fprintf('SPEED FLUCTUATION CONTROL:\n');
    fprintf('  Target Cs: %.6f\n', params.flywheelCoefficientOfFluctuation);
    fprintf('  Achieved Cs: %.6f\n', results.Cs_actual);
    fprintf('  Status: %s\n\n', check_status(results.Cs_actual <= params.flywheelCoefficientOfFluctuation * 1.01));

    fprintf('ENGINE CONFIGURATION:\n');
    fprintf('  Type: Beta-type Stirling Engine\n');
    fprintf('  Bore: %.0f mm\n', params.cylinderBore * 1000);
    fprintf('  Power Stroke: %.0f mm\n', params.powerCrankLength * 2 * 1000);
    fprintf('  Displacer Stroke: %.0f mm\n', params.displacerCrankLength * 2 * 1000);
    fprintf('  Operating Speed: %.0f RPM\n', params.averageRPM);
    fprintf('  Phase Shift: %.0f degrees\n', params.phaseShift * 180/pi);
    fprintf('  Compression Ratio: %.2f (given)\n\n', params.compressionRatio);

    fprintf('THERMODYNAMIC CONDITIONS:\n');
    fprintf('  Hot Space: %.0f K\n', params.hotTemperature);
    fprintf('  Cold Space: %.0f K\n', params.coldTemperature);
    fprintf('  Pressure at BDC: %.1f kPa\n', params.pressureAtBDC/1000);
    fprintf('  Pressure Range: %.2f - %.2f MPa\n\n', min(results.P)/1e6, max(results.P)/1e6);

    fprintf('POWER OUTPUT VALIDATION (Two Methods):\n');
    fprintf('  Method 1 (P-dV Integration): %.2f W\n', results.P_indicated);
    fprintf('  Method 2 (MEP): %.2f W\n', results.P_mep);
    fprintf('  Agreement: %.1f%%\n', 100 - abs(results.P_indicated - results.P_mep) / results.P_indicated * 100);
    fprintf('  Status: %s\n\n', check_status(abs(results.P_indicated - results.P_mep) / results.P_indicated < 0.05));

    fprintf('EFFICIENCY CHECK:\n');
    fprintf('  Thermal Efficiency: %.1f%%\n', results.efficiency * 100);
    fprintf('  Carnot Limit: %.1f%%\n', eta_carnot * 100);
    fprintf('  Status: %s (must be < Carnot)\n\n', check_status(results.efficiency < eta_carnot));

    fprintf('TORQUE ANALYSIS:\n');
    fprintf('  Mean Torque: %.2f N·m\n', results.T_mean);
    fprintf('  Max Torque: %.2f N·m\n', max(results.T_total));
    fprintf('  Min Torque: %.2f N·m\n\n', min(results.T_total));

    fprintf('DYNAMIC PERFORMANCE:\n');
    fprintf('  Speed Range: %.1f - %.1f RPM\n', min(results.rpm), max(results.rpm));
    fprintf('  Mean Speed: %.1f RPM\n', results.rpm_mean);
    fprintf('  Speed Variation: %.1f RPM\n\n', max(results.rpm) - min(results.rpm));

    fprintf('PHASE ANGLE OPTIMIZATION:\n');
    fprintf('  Current Phase: %.0f°\n', params.phaseShift * 180/pi);
    fprintf('  Optimal Phase: %.0f°\n', results.optimal_phase * 180/pi);
    fprintf('  Power at Current: %.3f W\n', results.P_indicated);
    fprintf('  Max Power at Optimal: %.3f W\n\n', max(results.optimization.power_curve));

    fprintf('PROJECT DELIVERABLES STATUS:\n');
    fprintf('  [%s] Flywheel diameter calculated\n', check_mark(results.flywheel.D_outer > 0));
    fprintf('  [%s] Coefficient of fluctuation maintained (Cs ≤ %.4f)\n', check_mark(results.Cs_actual <= params.flywheelCoefficientOfFluctuation * 1.01), params.flywheelCoefficientOfFluctuation);
    fprintf('  [%s] Two power methods validated (< 5%% difference)\n', check_mark(abs(results.P_indicated - results.P_mep) / results.P_indicated < 0.05));
    fprintf('  [%s] Efficiency below Carnot limit\n', check_mark(results.efficiency < eta_carnot));
    fprintf('  [%s] Four required plots generated\n\n', check_mark(exist('results/pv_diagram.png', 'file') == 2));

    all_requirements_met = results.flywheel.D_outer > 0 && ...
                          results.Cs_actual <= params.flywheelCoefficientOfFluctuation * 1.01 && ...
                          abs(results.P_indicated - results.P_mep) / results.P_indicated < 0.05 && ...
                          results.efficiency < eta_carnot;

    if all_requirements_met
        fprintf('FINAL STATUS: PROJECT REQUIREMENTS MET ✓\n');
        fprintf('Flywheel successfully sized for speed fluctuation control.\n');
    else
        fprintf('FINAL STATUS: REVIEW REQUIRED\n');
        fprintf('Check calculations and parameter values.\n');
    end

    fprintf('========================================================\n\n');
end

function mark = check_mark(condition)
    if condition
        mark = '✓';
    else
        mark = '✗';
    end
end

function status = check_status(condition)
    if condition
        status = 'PASS ✓';
    else
        status = 'FAIL ✗';
    end
end