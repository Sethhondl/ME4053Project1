function [optimal_phase, energy_curve, power_curve, efficiency_curve] = optimize_phase(params)
    original_phase = params.phaseShift;
    theta = linspace(0, 2*pi, params.simulationPointsPerCycle)';

    fprintf('\n=== PHASE ANGLE OPTIMIZATION (3 DECIMAL PRECISION) ===\n');

    %% Stage 1: Coarse Search (5° steps)
    fprintf('Stage 1: Coarse search (5° steps)...\n');
    coarse_range = 60:5:120;
    coarse_power = zeros(size(coarse_range));

    for i = 1:length(coarse_range)
        params_temp = params;
        params_temp.phaseShift = coarse_range(i) * pi/180;

        [V_total, V_exp, V_comp, ~, ~] = calc_volumes(theta, params_temp);
        [P, ~, ~] = schmidt_analysis(theta, V_total, V_exp, V_comp, params_temp);
        [~, P_indicated, ~, ~, ~, ~] = calc_power(P, V_total, theta, params_temp);

        coarse_power(i) = P_indicated;
    end

    [~, coarse_idx] = max(coarse_power);
    coarse_optimal = coarse_range(coarse_idx);
    fprintf('  Coarse optimum: %.0f°\n', coarse_optimal);

    %% Stage 2: Fine Search (0.1° steps)
    fprintf('Stage 2: Fine search (0.1° steps)...\n');
    fine_range = (coarse_optimal - 5):0.1:(coarse_optimal + 5);
    fine_power = zeros(size(fine_range));

    for i = 1:length(fine_range)
        params_temp = params;
        params_temp.phaseShift = fine_range(i) * pi/180;

        [V_total, V_exp, V_comp, ~, ~] = calc_volumes(theta, params_temp);
        [P, ~, ~] = schmidt_analysis(theta, V_total, V_exp, V_comp, params_temp);
        [~, P_indicated, ~, ~, ~, ~] = calc_power(P, V_total, theta, params_temp);

        fine_power(i) = P_indicated;
    end

    [~, fine_idx] = max(fine_power);
    fine_optimal = fine_range(fine_idx);
    fprintf('  Fine optimum: %.1f°\n', fine_optimal);

    %% Stage 3: Ultra-Fine Search (0.001° steps)
    fprintf('Stage 3: Ultra-fine search (0.001° steps)...\n');
    ultra_range = (fine_optimal - 0.1):0.001:(fine_optimal + 0.1);
    ultra_power = zeros(size(ultra_range));
    ultra_energy = zeros(size(ultra_range));
    ultra_efficiency = zeros(size(ultra_range));

    for i = 1:length(ultra_range)
        params_temp = params;
        params_temp.phaseShift = ultra_range(i) * pi/180;

        [V_total, V_exp, V_comp, ~, ~] = calc_volumes(theta, params_temp);
        [P, ~, ~] = schmidt_analysis(theta, V_total, V_exp, V_comp, params_temp);
        [W_indicated, P_indicated, ~, ~, ~, efficiency] = calc_power(P, V_total, theta, params_temp);

        ultra_power(i) = P_indicated;
        ultra_energy(i) = W_indicated;
        ultra_efficiency(i) = efficiency;
    end

    [max_power, ultra_idx] = max(ultra_power);
    ultra_optimal = ultra_range(ultra_idx);

    %% Stage 4: Parabolic Refinement for Maximum Precision
    fprintf('Stage 4: Parabolic refinement...\n');

    if ultra_idx > 1 && ultra_idx < length(ultra_range)
        x1 = ultra_range(ultra_idx - 1);
        x2 = ultra_range(ultra_idx);
        x3 = ultra_range(ultra_idx + 1);
        y1 = ultra_power(ultra_idx - 1);
        y2 = ultra_power(ultra_idx);
        y3 = ultra_power(ultra_idx + 1);

        denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
        A = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
        B = (x3^2 * (y1 - y2) + x2^2 * (y3 - y1) + x1^2 * (y2 - y3)) / denom;

        if A < 0
            parabolic_optimal = -B / (2 * A);

            params_temp = params;
            params_temp.phaseShift = parabolic_optimal * pi/180;
            [V_total, V_exp, V_comp, ~, ~] = calc_volumes(theta, params_temp);
            [P, ~, ~] = schmidt_analysis(theta, V_total, V_exp, V_comp, params_temp);
            [~, P_parabolic, ~, ~, ~, ~] = calc_power(P, V_total, theta, params_temp);

            if P_parabolic > max_power
                optimal_phase_deg = parabolic_optimal;
                max_power = P_parabolic;
                fprintf('  Parabolic refinement: %.3f° (improved)\n', optimal_phase_deg);
            else
                optimal_phase_deg = ultra_optimal;
                fprintf('  Parabolic refinement: %.3f° (no improvement)\n', ultra_optimal);
            end
        else
            optimal_phase_deg = ultra_optimal;
            fprintf('  Parabolic refinement skipped (invalid curvature)\n');
        end
    else
        optimal_phase_deg = ultra_optimal;
    end

    %% Compile results for all evaluation points
    all_phases = [coarse_range, fine_range, ultra_range];
    all_power = [coarse_power, fine_power, ultra_power];

    all_energy = zeros(size(all_phases));
    all_efficiency = zeros(size(all_phases));

    unique_phases = unique(all_phases);
    power_curve = zeros(size(unique_phases));
    energy_curve = zeros(size(unique_phases));
    efficiency_curve = zeros(size(unique_phases));

    for i = 1:length(unique_phases)
        params_temp = params;
        params_temp.phaseShift = unique_phases(i) * pi/180;

        [V_total, V_exp, V_comp, ~, ~] = calc_volumes(theta, params_temp);
        [P, ~, ~] = schmidt_analysis(theta, V_total, V_exp, V_comp, params_temp);
        [W_indicated, P_indicated, ~, ~, ~, efficiency] = calc_power(P, V_total, theta, params_temp);

        power_curve(i) = P_indicated;
        energy_curve(i) = W_indicated;
        efficiency_curve(i) = efficiency;
    end

    %% Final validation at optimal point
    params_temp = params;
    params_temp.phaseShift = optimal_phase_deg * pi/180;
    [V_total, V_exp, V_comp, ~, ~] = calc_volumes(theta, params_temp);
    [P, ~, ~] = schmidt_analysis(theta, V_total, V_exp, V_comp, params_temp);
    [W_opt, P_opt, ~, ~, ~, eff_opt] = calc_power(P, V_total, theta, params_temp);

    fprintf('\n=== OPTIMIZATION RESULTS ===\n');
    fprintf('Optimal phase angle: %.3f°\n', optimal_phase_deg);
    fprintf('Power at optimum: %.3f W\n', P_opt);
    fprintf('Efficiency at optimum: %.3f%%\n', eff_opt * 100);
    fprintf('Energy per cycle: %.3f J\n', W_opt);

    original_phase_deg = original_phase * 180/pi;
    [~, orig_idx] = min(abs(unique_phases - original_phase_deg));
    if orig_idx <= length(unique_phases)
        improvement = (P_opt - power_curve(orig_idx)) / power_curve(orig_idx) * 100;
        fprintf('\nImprovement over %.0f°: %.2f%%\n', original_phase_deg, improvement);
    end

    optimal_phase = optimal_phase_deg;

    params.phaseShift = original_phase;
end