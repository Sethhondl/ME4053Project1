function [optimal_phase, energy_curve, power_curve, efficiency_curve] = optimize_phase(params)
    % OPTIMIZE_PHASE Find optimal phase angle using adaptive step size
    % Uses multi-stage optimization: coarse → medium → fine → ultra-fine
    %
    % Inputs:
    %   params - engine parameters structure
    %
    % Outputs:
    %   optimal_phase - optimal phase angle in degrees
    %   energy_curve - energy per cycle vs phase angle [phase, energy]
    %   power_curve - power output vs phase angle [phase, power]
    %   efficiency_curve - thermal efficiency vs phase angle [phase, efficiency]

    % Save original phase angle
    original_phase = params.phaseShift;

    % Simulation parameters
    theta = linspace(0, 2*pi, params.simulationPointsPerCycle);

    fprintf('\n========== ADAPTIVE PHASE ANGLE OPTIMIZATION ==========\n');
    fprintf('Using multi-stage search with adaptive step sizes\n\n');

    % Initialize storage for all evaluations
    all_phases = [];
    all_power = [];
    all_energy = [];
    all_efficiency = [];
    all_torque_rms = [];
    all_pressure_ratio = [];
    all_cs = [];

    %% STAGE 1: COARSE SEARCH (10° steps)
    fprintf('STAGE 1: Coarse Search (10° steps)\n');
    fprintf('----------------------------------------\n');
    coarse_range = 30:10:150;  % Wider initial range
    n_coarse = length(coarse_range);

    coarse_power = zeros(n_coarse, 1);
    coarse_energy = zeros(n_coarse, 1);
    coarse_efficiency = zeros(n_coarse, 1);

    for i = 1:n_coarse
        phase = coarse_range(i);
        params.phaseShift = phase * pi/180;

        % Calculate performance
        [V_total, V_exp, V_comp, x_power, x_disp] = calc_volumes(theta, params);
        [P, ~, P_mean] = schmidt_analysis(theta, V_total, V_exp, V_comp, params);
        [W_indicated, P_indicated, ~, ~, ~, efficiency] = calc_power(P, V_total, theta, params);

        coarse_power(i) = P_indicated;
        coarse_energy(i) = W_indicated;
        coarse_efficiency(i) = efficiency;

        if mod(phase, 30) == 0  % Print every 30 degrees
            fprintf('  φ = %3d°: P = %6.2f W, η = %4.1f%%\n', ...
                    phase, P_indicated, efficiency*100);
        end
    end

    % Find coarse optimum
    [max_coarse_power, coarse_idx] = max(coarse_power);
    coarse_optimal = coarse_range(coarse_idx);

    fprintf('  → Coarse optimum: %d° (P = %.2f W)\n\n', ...
            coarse_optimal, max_coarse_power);

    % Store coarse results
    all_phases = [all_phases, coarse_range];
    all_power = [all_power, coarse_power'];
    all_energy = [all_energy, coarse_energy'];
    all_efficiency = [all_efficiency, coarse_efficiency'];

    %% STAGE 2: MEDIUM SEARCH (2° steps around coarse optimum)
    fprintf('STAGE 2: Medium Search (2° steps)\n');
    fprintf('----------------------------------------\n');
    medium_center = coarse_optimal;
    medium_range = max(30, medium_center-20):2:min(150, medium_center+20);
    n_medium = length(medium_range);

    medium_power = zeros(n_medium, 1);
    medium_energy = zeros(n_medium, 1);
    medium_efficiency = zeros(n_medium, 1);
    medium_torque_rms = zeros(n_medium, 1);
    medium_pressure_ratio = zeros(n_medium, 1);

    for i = 1:n_medium
        phase = medium_range(i);
        params.phaseShift = phase * pi/180;

        % Calculate performance with additional metrics
        [V_total, V_exp, V_comp, x_power, x_disp] = calc_volumes(theta, params);
        [P, ~, P_mean] = schmidt_analysis(theta, V_total, V_exp, V_comp, params);
        [W_indicated, P_indicated, ~, ~, ~, efficiency] = calc_power(P, V_total, theta, params);

        % Additional metrics
        [T_total, ~, ~, T_mean] = calc_torque(P, theta, x_power, x_disp, params);
        T_rms = sqrt(mean(T_total.^2));
        P_ratio = max(P)/min(P);

        medium_power(i) = P_indicated;
        medium_energy(i) = W_indicated;
        medium_efficiency(i) = efficiency;
        medium_torque_rms(i) = T_rms;
        medium_pressure_ratio(i) = P_ratio;

        if mod(phase, 10) == 0  % Print every 10 degrees
            fprintf('  φ = %3d°: P = %6.2f W, η = %4.1f%%, T_rms = %5.2f Nm\n', ...
                    phase, P_indicated, efficiency*100, T_rms);
        end
    end

    % Find medium optimum
    [max_medium_power, medium_idx] = max(medium_power);
    medium_optimal = medium_range(medium_idx);

    fprintf('  → Medium optimum: %d° (P = %.2f W)\n\n', ...
            medium_optimal, max_medium_power);

    % Store medium results (avoid duplicates)
    new_phases = setdiff(medium_range, all_phases);
    if ~isempty(new_phases)
        [~, new_indices] = ismember(new_phases, medium_range);
        all_phases = [all_phases, new_phases];
        all_power = [all_power, medium_power(new_indices)'];
        all_energy = [all_energy, medium_energy(new_indices)'];
        all_efficiency = [all_efficiency, medium_efficiency(new_indices)'];
        all_torque_rms = [all_torque_rms, medium_torque_rms(new_indices)'];
        all_pressure_ratio = [all_pressure_ratio, medium_pressure_ratio(new_indices)'];
    end

    %% STAGE 3: FINE SEARCH (0.5° steps around medium optimum)
    fprintf('STAGE 3: Fine Search (0.5° steps)\n');
    fprintf('----------------------------------------\n');
    fine_center = medium_optimal;
    fine_range = max(30, fine_center-5):0.5:min(150, fine_center+5);
    n_fine = length(fine_range);

    fine_power = zeros(n_fine, 1);
    fine_energy = zeros(n_fine, 1);
    fine_efficiency = zeros(n_fine, 1);
    fine_torque_rms = zeros(n_fine, 1);
    fine_pressure_ratio = zeros(n_fine, 1);
    fine_cs = zeros(n_fine, 1);

    for i = 1:n_fine
        phase = fine_range(i);
        params.phaseShift = phase * pi/180;

        % Calculate comprehensive performance metrics
        [V_total, V_exp, V_comp, x_power, x_disp] = calc_volumes(theta, params);
        [P, ~, P_mean] = schmidt_analysis(theta, V_total, V_exp, V_comp, params);
        [W_indicated, P_indicated, ~, ~, MEP, efficiency] = calc_power(P, V_total, theta, params);

        % Torque analysis
        [T_total, T_power, T_disp, T_mean] = calc_torque(P, theta, x_power, x_disp, params);
        T_rms = sqrt(mean(T_total.^2));

        % Pressure metrics
        P_ratio = max(P)/min(P);

        % Speed fluctuation (simplified calculation)
        if T_mean > 0
            energy_var = max(cumtrapz(theta, T_total)) - min(cumtrapz(theta, T_total));
            I_est = energy_var / (params.averageAngularVelocity^2 * params.flywheelCoefficientOfFluctuation);
            [omega, ~, ~, Cs] = simulate_dynamics(T_total, theta, I_est, params);
            fine_cs(i) = Cs;
        else
            fine_cs(i) = NaN;
        end

        fine_power(i) = P_indicated;
        fine_energy(i) = W_indicated;
        fine_efficiency(i) = efficiency;
        fine_torque_rms(i) = T_rms;
        fine_pressure_ratio(i) = P_ratio;

        if mod(phase, 2) == 0  % Print every 2 degrees
            fprintf('  φ = %5.1f°: P = %6.2f W, η = %4.1f%%, MEP = %5.1f kPa, Cs = %.4f\n', ...
                    phase, P_indicated, efficiency*100, MEP/1000, fine_cs(i));
        end
    end

    % Find fine optimum
    [max_fine_power, fine_idx] = max(fine_power);
    fine_optimal = fine_range(fine_idx);

    fprintf('  → Fine optimum: %.1f° (P = %.2f W)\n\n', ...
            fine_optimal, max_fine_power);

    %% STAGE 4: ULTRA-FINE SEARCH (0.01° steps for better precision)
    fprintf('STAGE 4: Ultra-Fine Search (0.01° steps)\n');
    fprintf('----------------------------------------\n');
    ultra_center = fine_optimal;
    ultra_range = max(30, ultra_center-0.5):0.01:min(150, ultra_center+0.5);
    n_ultra = length(ultra_range);

    ultra_power = zeros(n_ultra, 1);
    ultra_energy = zeros(n_ultra, 1);
    ultra_efficiency = zeros(n_ultra, 1);
    ultra_mep = zeros(n_ultra, 1);

    best_power = 0;
    best_phase = ultra_center;

    for i = 1:n_ultra
        phase = ultra_range(i);
        params.phaseShift = phase * pi/180;

        % Precise calculation
        [V_total, V_exp, V_comp, ~, ~] = calc_volumes(theta, params);
        [P, ~, ~] = schmidt_analysis(theta, V_total, V_exp, V_comp, params);
        [W_indicated, P_indicated, ~, ~, MEP, efficiency] = calc_power(P, V_total, theta, params);

        ultra_power(i) = P_indicated;
        ultra_energy(i) = W_indicated;
        ultra_efficiency(i) = efficiency;
        ultra_mep(i) = MEP;

        if P_indicated > best_power
            best_power = P_indicated;
            best_phase = phase;
        end
    end

    % Ultra-fine optimum
    [max_ultra_power, ultra_idx] = max(ultra_power);
    ultra_optimal = ultra_range(ultra_idx);

    fprintf('  → Ultra-fine optimum: %.2f° (P = %.2f W)\n\n', ...
            ultra_optimal, max_ultra_power);

    %% STAGE 5: PRECISION SEARCH (0.001° steps for 3 decimal precision)
    fprintf('STAGE 5: Precision Search (0.001° steps)\n');
    fprintf('----------------------------------------\n');
    precision_center = ultra_optimal;
    precision_range = max(30, precision_center-0.05):0.001:min(150, precision_center+0.05);
    n_precision = length(precision_range);

    precision_power = zeros(n_precision, 1);
    precision_energy = zeros(n_precision, 1);
    precision_efficiency = zeros(n_precision, 1);
    precision_mep = zeros(n_precision, 1);

    for i = 1:n_precision
        phase = precision_range(i);
        params.phaseShift = phase * pi/180;

        % High precision calculation
        [V_total, V_exp, V_comp, ~, ~] = calc_volumes(theta, params);
        [P, ~, ~] = schmidt_analysis(theta, V_total, V_exp, V_comp, params);
        [W_indicated, P_indicated, ~, ~, MEP, efficiency] = calc_power(P, V_total, theta, params);

        precision_power(i) = P_indicated;
        precision_energy(i) = W_indicated;
        precision_efficiency(i) = efficiency;
        precision_mep(i) = MEP;
    end

    % Find precision optimum
    [max_precision_power, precision_idx] = max(precision_power);
    precision_optimal = precision_range(precision_idx);

    %% STAGE 6: Parabolic Interpolation for Final Refinement
    fprintf('STAGE 6: Parabolic Interpolation\n');
    fprintf('----------------------------------------\n');

    % Use parabolic interpolation around the peak for sub-grid precision
    if precision_idx > 1 && precision_idx < n_precision
        % Get three points around maximum
        x1 = precision_range(precision_idx - 1);
        x2 = precision_range(precision_idx);
        x3 = precision_range(precision_idx + 1);
        y1 = precision_power(precision_idx - 1);
        y2 = precision_power(precision_idx);
        y3 = precision_power(precision_idx + 1);

        % Fit parabola through three points
        denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
        A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
        B = (x3^2 * (y1 - y2) + x2^2 * (y3 - y1) + x1^2 * (y2 - y3)) / denom;

        % Find vertex of parabola (maximum if A < 0)
        if A < 0
            parabolic_optimal = -B / (2 * A);

            % Verify this is actually better
            params.phaseShift = parabolic_optimal * pi/180;
            [V_total, V_exp, V_comp, ~, ~] = calc_volumes(theta, params);
            [P, ~, ~] = schmidt_analysis(theta, V_total, V_exp, V_comp, params);
            [~, P_parabolic, ~, ~, ~, ~] = calc_power(P, V_total, theta, params);

            if P_parabolic > max_precision_power
                optimal_phase = parabolic_optimal;
                max_final_power = P_parabolic;
                fprintf('  → Parabolic optimum: %.3f° (P = %.3f W) [IMPROVED]\n\n', ...
                        optimal_phase, max_final_power);
            else
                optimal_phase = precision_optimal;
                max_final_power = max_precision_power;
                fprintf('  → Parabolic optimum: %.3f° (P = %.3f W) [NO IMPROVEMENT]\n\n', ...
                        precision_optimal, max_precision_power);
            end
        else
            optimal_phase = precision_optimal;
            max_final_power = max_precision_power;
            fprintf('  → Grid optimum retained: %.3f° (invalid parabola)\n\n', ...
                    precision_optimal, max_precision_power);
        end
    else
        optimal_phase = precision_optimal;
        max_final_power = max_precision_power;
        fprintf('  → Final optimum: %.3f° (P = %.3f W)\n\n', ...
                precision_optimal, max_precision_power);
    end

    %% COMPILE ALL RESULTS
    % Sort all results by phase angle
    [all_phases_sorted, sort_idx] = sort(all_phases);
    all_power_sorted = all_power(sort_idx);
    all_energy_sorted = all_energy(sort_idx);
    all_efficiency_sorted = all_efficiency(sort_idx);

    % Add ultra-fine and precision results
    all_phases_final = [all_phases_sorted, ultra_range, precision_range];
    all_power_final = [all_power_sorted, ultra_power', precision_power'];
    all_energy_final = [all_energy_sorted, ultra_energy', precision_energy'];
    all_efficiency_final = [all_efficiency_sorted, ultra_efficiency', precision_efficiency'];

    % Remove duplicates and sort
    [unique_phases, unique_idx] = unique(all_phases_final);
    power_curve = all_power_final(unique_idx);
    energy_curve = all_energy_final(unique_idx);
    efficiency_curve = all_efficiency_final(unique_idx);

    %% ANALYSIS SUMMARY
    fprintf('========== OPTIMIZATION SUMMARY ==========\n');
    fprintf('Total evaluations: %d\n', length(unique_phases));
    fprintf('  Stage 1 (coarse): %d evaluations\n', n_coarse);
    fprintf('  Stage 2 (medium): %d evaluations\n', n_medium);
    fprintf('  Stage 3 (fine): %d evaluations\n', n_fine);
    fprintf('  Stage 4 (ultra): %d evaluations\n', n_ultra);
    fprintf('  Stage 5 (precision): %d evaluations\n', n_precision);
    fprintf('\nConvergence:\n');
    fprintf('  Coarse → Medium: %.0f° → %.0f°\n', coarse_optimal, medium_optimal);
    fprintf('  Medium → Fine: %.0f° → %.1f°\n', medium_optimal, fine_optimal);
    fprintf('  Fine → Ultra: %.1f° → %.2f°\n', fine_optimal, ultra_optimal);
    fprintf('  Ultra → Precision: %.2f° → %.3f°\n', ultra_optimal, optimal_phase);

    % Find maximum power and efficiency
    [max_power, max_p_idx] = max(power_curve);
    [max_eff, max_e_idx] = max(efficiency_curve);
    phase_max_power = unique_phases(max_p_idx);
    phase_max_eff = unique_phases(max_e_idx);

    fprintf('\n========== FINAL RESULTS ==========\n');
    fprintf('Optimal Phase (Max Power): %.3f°\n', optimal_phase);
    fprintf('  Power: %.3f W (%.3f kW)\n', max_final_power, max_final_power/1000);
    fprintf('  Efficiency: %.3f%%\n', efficiency_curve(max_p_idx)*100);
    fprintf('  Energy/cycle: %.3f J\n', energy_curve(max_p_idx));

    fprintf('\nOptimal Phase (Max Efficiency): %.3f°\n', phase_max_eff);
    fprintf('  Efficiency: %.3f%%\n', max_eff*100);
    fprintf('  Power: %.3f W (%.3f kW)\n', power_curve(max_e_idx), power_curve(max_e_idx)/1000);

    % Compare with original setting
    original_phase_deg = original_phase * 180/pi;
    [~, orig_idx] = min(abs(unique_phases - original_phase_deg));
    if orig_idx <= length(unique_phases)
        orig_power = power_curve(orig_idx);
        orig_eff = efficiency_curve(orig_idx);
        power_improvement = (max_power - orig_power) / orig_power * 100;

        fprintf('\nComparison with Original (%.0f°):\n', original_phase_deg);
        fprintf('  Original Power: %.2f W\n', orig_power);
        fprintf('  Original Efficiency: %.1f%%\n', orig_eff*100);
        if power_improvement > 0
            fprintf('  Potential Power Gain: %.1f%%\n', power_improvement);
        else
            fprintf('  Current setting is near optimal\n');
        end
    end

    % Restore original phase angle
    params.phaseShift = original_phase;

    % Return results as matrices with phase angle in first column
    energy_curve = [unique_phases(:), energy_curve(:)];
    power_curve = [unique_phases(:), power_curve(:)];
    efficiency_curve = [unique_phases(:), efficiency_curve(:)];

    fprintf('=====================================\n\n');
end