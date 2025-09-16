function [optimal_phase, energy_curve, power_curve, efficiency_curve] = optimize_phase(params)
    % OPTIMIZE_PHASE Find optimal phase angle for maximum energy/power
    %
    % Inputs:
    %   params - engine parameters structure
    %
    % Outputs:
    %   optimal_phase - optimal phase angle in degrees
    %   energy_curve - energy per cycle vs phase angle
    %   power_curve - power output vs phase angle
    %   efficiency_curve - thermal efficiency vs phase angle
    
    % Define phase angle range to test (45 to 135 degrees)
    % Literature shows optimal phases can be 57.6°, 80.78°, 90°, 103° etc.
    % Expanding range to capture all possible optimal values
    phase_range = 45:5:135;  % degrees - expanded range based on research
    n_phases = length(phase_range);
    
    % Initialize arrays for results
    energy_curve = zeros(n_phases, 1);
    power_curve = zeros(n_phases, 1);
    efficiency_curve = zeros(n_phases, 1);
    
    % Save original phase angle
    original_phase = params.phase_shift;
    
    % Simulation parameters
    theta = linspace(0, 2*pi, params.sim.n_points);
    
    fprintf('Phase Angle Optimization:\n');
    fprintf('Testing phase angles from %d to %d degrees...\n', ...
            min(phase_range), max(phase_range));
    
    % Test each phase angle
    for i = 1:n_phases
        % Update phase angle
        params.phase_shift = phase_range(i) * pi/180;  % Convert to radians
        
        % Calculate volumes
        [V_total, V_exp, V_comp, ~, ~] = calc_volumes(theta, params);
        
        % Calculate pressure (Schmidt analysis)
        [P, ~, ~] = schmidt_analysis(theta, V_total, V_exp, V_comp, params);
        
        % Calculate work and power
        [W_indicated, P_indicated, ~, ~, ~, efficiency] = ...
            calc_power(P, V_total, theta, params);
        
        % Store results
        energy_curve(i) = W_indicated;
        power_curve(i) = P_indicated;
        efficiency_curve(i) = efficiency;
        
        % Progress indicator
        if mod(i, 5) == 0
            fprintf('  Phase = %d deg: Power = %.2f kW, Efficiency = %.1f%%\n', ...
                    phase_range(i), P_indicated/1000, efficiency*100);
        end
    end
    
    % Find optimal phase angle (based on power)
    [max_power, max_idx] = max(power_curve);
    optimal_phase = phase_range(max_idx);
    optimal_energy = energy_curve(max_idx);
    optimal_efficiency = efficiency_curve(max_idx);
    
    % Also find phase for maximum efficiency
    [max_eff, eff_idx] = max(efficiency_curve);
    phase_max_eff = phase_range(eff_idx);
    
    % Restore original phase angle
    params.phase_shift = original_phase;
    
    % Display optimization results
    fprintf('\nOptimization Results:\n');
    fprintf('  Optimal Phase Angle (Max Power): %d degrees\n', optimal_phase);
    fprintf('    Power at Optimum: %.2f kW\n', max_power/1000);
    fprintf('    Energy per Cycle: %.2f J\n', optimal_energy);
    fprintf('    Efficiency: %.1f%%\n', optimal_efficiency*100);
    fprintf('  Phase for Max Efficiency: %d degrees\n', phase_max_eff);
    fprintf('    Max Efficiency: %.1f%%\n', max_eff*100);
    fprintf('    Power at Max Efficiency: %.2f kW\n', ...
            power_curve(eff_idx)/1000);
    
    % Compare with current setting
    current_phase_deg = original_phase * 180/pi;
    [~, current_idx] = min(abs(phase_range - current_phase_deg));
    if current_idx <= n_phases
        current_power = power_curve(current_idx);
        power_improvement = (max_power - current_power) / current_power * 100;
        fprintf('  Current Phase: %.0f degrees\n', current_phase_deg);
        fprintf('    Current Power: %.2f kW\n', current_power/1000);
        if power_improvement > 0
            fprintf('    Potential Improvement: %.1f%%\n', power_improvement);
        end
    end
    
    fprintf('\n');
    
    % Return phase angles in degrees for plotting
    energy_curve = [phase_range(:), energy_curve];
    power_curve = [phase_range(:), power_curve];
    efficiency_curve = [phase_range(:), efficiency_curve];
end