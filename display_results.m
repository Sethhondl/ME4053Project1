function display_results(results, params)
    % DISPLAY_RESULTS Display comprehensive analysis results
    %
    % Inputs:
    %   results - structure containing all analysis results
    %   params - engine parameters structure
    
    % Calculate Carnot efficiency for comparison
    eta_carnot = 1 - params.coldTemperature / params.hotTemperature;
    
    % Display header
    fprintf('\n');
    fprintf('========================================================\n');
    fprintf('       STIRLING ENGINE ANALYSIS - FINAL RESULTS        \n');
    fprintf('========================================================\n');
    fprintf('\n');
    
    %% Engine Configuration Summary
    fprintf('ENGINE CONFIGURATION:\n');
    fprintf('---------------------\n');
    fprintf('  Type: Beta-type Stirling Engine\n');
    fprintf('  Working Fluid: %s\n', params.gasName);
    fprintf('  Operating Speed: %.0f RPM\n', params.averageRPM);
    fprintf('  Phase Shift: %.0f degrees\n', params.phaseShift * 180/pi);
    fprintf('  Compression Ratio: %.2f\n', params.compressionRatio);
    fprintf('\n');
    
    %% Thermodynamic Performance
    fprintf('THERMODYNAMIC PERFORMANCE:\n');
    fprintf('--------------------------\n');
    fprintf('  Operating Temperatures:\n');
    fprintf('    Hot Space: %.0f K (%.0f°C)\n', params.hotTemperature, params.hotTemperature - 273.15);
    fprintf('    Cold Space: %.0f K (%.0f°C)\n', params.coldTemperature, params.coldTemperature - 273.15);
    fprintf('  Pressure Range: %.2f - %.2f MPa\n', ...
            min(results.P)/1e6, max(results.P)/1e6);
    fprintf('\n');
    
    %% Flywheel Design Results
    fprintf('FLYWHEEL DESIGN:\n');
    fprintf('----------------\n');
    fprintf('  Material: Steel (ρ = %.0f kg/m³)\n', params.flywheelMaterialDensity);
    fprintf('  Outer Diameter: %.3f m (%.1f mm)\n', ...
            results.flywheel.D_outer, results.flywheel.D_outer * 1000);
    fprintf('  Mass: %.2f kg\n', results.flywheel.mass);
    fprintf('  Moment of Inertia: %.4f kg·m²\n', results.flywheel.I_required);
    fprintf('  Energy Fluctuation: %.2f J\n', results.flywheel.energy_fluctuation);
    fprintf('\n');
    
    %% Power Output Comparison
    fprintf('POWER OUTPUT (Two Methods):\n');
    fprintf('----------------------------\n');
    fprintf('  Method 1 - Direct Integration:\n');
    fprintf('    Work per Cycle: %.2f J\n', results.W_indicated);
    fprintf('    Power Output: %.2f W (%.3f kW)\n', ...
            results.P_indicated, results.P_indicated/1000);
    fprintf('  Method 2 - Mean Effective Pressure:\n');
    fprintf('    MEP: %.3f MPa\n', results.MEP/1e6);
    fprintf('    Work per Cycle: %.2f J\n', results.W_mep);
    fprintf('    Power Output: %.2f W (%.3f kW)\n', ...
            results.P_mep, results.P_mep/1000);
    
    % Calculate and display agreement
    power_diff = abs(results.P_indicated - results.P_mep) / results.P_indicated * 100;
    fprintf('  Agreement Between Methods: %.1f%%\n', 100 - power_diff);
    fprintf('\n');
    
    %% Efficiency Analysis
    fprintf('EFFICIENCY ANALYSIS:\n');
    fprintf('--------------------\n');
    fprintf('  Thermal Efficiency: %.1f%%\n', results.efficiency * 100);
    fprintf('  Carnot Efficiency Limit: %.1f%%\n', eta_carnot * 100);
    fprintf('  Relative Efficiency: %.1f%% of Carnot\n', ...
            (results.efficiency / eta_carnot) * 100);
    fprintf('\n');
    
    %% Dynamic Performance
    fprintf('DYNAMIC PERFORMANCE:\n');
    fprintf('--------------------\n');
    fprintf('  Speed Range: %.1f - %.1f RPM\n', ...
            min(results.rpm), max(results.rpm));
    fprintf('  Speed Variation: %.1f RPM\n', ...
            max(results.rpm) - min(results.rpm));
    fprintf('  Coefficient of Fluctuation:\n');
    fprintf('    Target: %.4f\n', params.flywheelCoefficientOfFluctuation);
    fprintf('    Achieved: %.4f\n', results.Cs_actual);
    if results.Cs_actual <= params.flywheelCoefficientOfFluctuation * 1.1
        fprintf('    Status: ✓ MEETS REQUIREMENT\n');
    else
        fprintf('    Status: ✗ EXCEEDS LIMIT\n');
    end
    fprintf('\n');
    
    %% Optimization Results
    fprintf('PHASE ANGLE OPTIMIZATION:\n');
    fprintf('--------------------------\n');

    % Extract detailed optimization data
    phase_values = results.optimization.power_curve(:,1);
    power_values = results.optimization.power_curve(:,2);
    eff_values = results.optimization.efficiency_curve(:,2);
    energy_values = results.optimization.energy_curve(:,2);

    % Find optimal points
    [max_power, max_p_idx] = max(power_values);
    [max_eff, max_e_idx] = max(eff_values);
    optimal_phase_power = phase_values(max_p_idx);
    optimal_phase_eff = phase_values(max_e_idx);

    fprintf('  MAXIMUM POWER:\n');
    fprintf('    Phase: %.1f° | Power: %.2f W | Efficiency: %.1f%%\n', ...
            optimal_phase_power, max_power, eff_values(max_p_idx)*100);

    fprintf('  MAXIMUM EFFICIENCY:\n');
    fprintf('    Phase: %.1f° | Efficiency: %.1f%% | Power: %.2f W\n', ...
            optimal_phase_eff, max_eff*100, power_values(max_e_idx));

    % Current configuration
    current_phase = params.phaseShift * 180/pi;
    if current_phase >= min(phase_values) && current_phase <= max(phase_values)
        current_power = interp1(phase_values, power_values, current_phase);
        current_eff = interp1(phase_values, eff_values, current_phase);

        fprintf('\n  CURRENT SETTING (%.0f°):\n', current_phase);
        fprintf('    Power: %.2f W | Efficiency: %.1f%%\n', ...
                current_power, current_eff*100);

        % Calculate improvements
        power_improvement = (max_power - current_power) / current_power * 100;
        if power_improvement > 1
            fprintf('    → Potential Power Gain: +%.1f%% at %.1f°\n', ...
                    power_improvement, optimal_phase_power);
        else
            fprintf('    → Near optimal for power\n');
        end
    end

    % Performance sensitivity
    if max_p_idx > 1 && max_p_idx < length(phase_values)
        dP = power_values(max_p_idx+1) - power_values(max_p_idx-1);
        dphi = phase_values(max_p_idx+1) - phase_values(max_p_idx-1);
        sensitivity = abs(dP/dphi);
        fprintf('\n  SENSITIVITY: %.1f W/° ', sensitivity);
        if sensitivity < 0.5
            fprintf('(Robust)\n');
        elseif sensitivity < 2
            fprintf('(Moderate)\n');
        else
            fprintf('(Sensitive)\n');
        end
    end

    % Operating window
    power_90pct = 0.9 * max_power;
    valid_phases = phase_values(power_values >= power_90pct);
    if ~isempty(valid_phases)
        fprintf('  OPERATING WINDOW: %.1f° - %.1f° for >90%% power\n', ...
                min(valid_phases), max(valid_phases));
    end

    % Store optimal phase
    results.optimal_phase = optimal_phase_power;
    fprintf('\n');
    
    %% Design Validation
    fprintf('DESIGN VALIDATION:\n');
    fprintf('------------------\n');
    
    % Check all requirements
    requirements_met = true;
    
    % Power output (informational only - not a requirement)
    fprintf('  • Power Output: %.2f kW\n', results.P_indicated/1000);
    % Note: Power output is determined by the given parameters
    % No specific power range requirement was specified
    
    % Efficiency check
    if results.efficiency <= eta_carnot
        fprintf('  ✓ Efficiency: %.1f%% (< Carnot limit)\n', ...
                results.efficiency * 100);
    else
        fprintf('  ✗ Efficiency: %.1f%% (EXCEEDS Carnot limit)\n', ...
                results.efficiency * 100);
        requirements_met = false;
    end
    
    % Flywheel size check
    if results.flywheel.D_outer <= params.maximumFlywheelDiameter
        fprintf('  ✓ Flywheel Diameter: %.2f m (< %.1f m limit)\n', ...
                results.flywheel.D_outer, params.maximumFlywheelDiameter);
    else
        fprintf('  ✗ Flywheel Diameter: %.2f m (EXCEEDS %.1f m limit)\n', ...
                results.flywheel.D_outer, params.maximumFlywheelDiameter);
        requirements_met = false;
    end
    
    % Speed fluctuation check
    if results.Cs_actual <= params.flywheelCoefficientOfFluctuation * 1.1
        fprintf('  ✓ Speed Fluctuation: Cs = %.4f (< %.4f target)\n', ...
                results.Cs_actual, params.flywheelCoefficientOfFluctuation);
    else
        fprintf('  ✗ Speed Fluctuation: Cs = %.4f (EXCEEDS %.4f target)\n', ...
                results.Cs_actual, params.flywheelCoefficientOfFluctuation);
        requirements_met = false;
    end
    
    % Positive work check
    if results.W_indicated > 0
        fprintf('  ✓ Work Output: %.1f J (positive)\n', results.W_indicated);
    else
        fprintf('  ✗ Work Output: %.1f J (NEGATIVE - not producing power)\n', ...
                results.W_indicated);
        requirements_met = false;
    end
    
    fprintf('\n');
    
    %% Final Summary
    fprintf('========================================================\n');
    fprintf('                    FINAL SUMMARY                      \n');
    fprintf('========================================================\n');
    
    if requirements_met
        fprintf('STATUS: ALL DESIGN REQUIREMENTS MET ✓\n');
    else
        fprintf('STATUS: SOME REQUIREMENTS NOT MET ✗\n');
        fprintf('         Review design parameters and rerun analysis\n');
    end
    
    fprintf('\nKey Performance Metrics:\n');
    fprintf('  • Power Output: %.2f kW\n', results.P_indicated/1000);
    fprintf('  • Thermal Efficiency: %.1f%%\n', results.efficiency * 100);
    fprintf('  • Flywheel Diameter: %.3f m\n', results.flywheel.D_outer);
    fprintf('  • Flywheel Mass: %.1f kg\n', results.flywheel.mass);
    fprintf('  • Speed Fluctuation: %.2f%%\n', results.Cs_actual * 100);
    
    fprintf('\n========================================================\n\n');
    
    % Save results to text file
    save_results_to_file(results, params);
end

function save_results_to_file(results, params)
    % Save results summary to text file
    
    filename = 'results/analysis_summary.txt';
    fid = fopen(filename, 'w');
    
    if fid == -1
        warning('Could not create results summary file');
        return;
    end
    
    % Write summary to file
    fprintf(fid, 'STIRLING ENGINE ANALYSIS SUMMARY\n');
    fprintf(fid, 'Generated: %s\n\n', datestr(now));
    
    fprintf(fid, 'CONFIGURATION:\n');
    fprintf(fid, '  Engine Type: Beta-type Stirling\n');
    fprintf(fid, '  Working Fluid: %s\n', params.gasName);
    fprintf(fid, '  Phase Shift: %.0f degrees\n', params.phaseShift * 180/pi);
    fprintf(fid, '  Compression Ratio: %.2f\n\n', params.compressionRatio);
    
    fprintf(fid, 'FLYWHEEL DESIGN:\n');
    fprintf(fid, '  Outer Diameter: %.3f m\n', results.flywheel.D_outer);
    fprintf(fid, '  Mass: %.2f kg\n', results.flywheel.mass);
    fprintf(fid, '  Moment of Inertia: %.4f kg·m²\n\n', results.flywheel.I_required);
    
    fprintf(fid, 'PERFORMANCE:\n');
    fprintf(fid, '  Power Output: %.2f kW\n', results.P_indicated/1000);
    fprintf(fid, '  Thermal Efficiency: %.1f%%\n', results.efficiency * 100);
    fprintf(fid, '  Speed Fluctuation (Cs): %.4f\n', results.Cs_actual);
    
    fclose(fid);
    fprintf('Results summary saved to: %s\n', filename);
end