function [optimal_phase, energy_curve, power_curve, efficiency_curve] = optimize_phase(params)
    % Find optimal phase angle for maximum power

    phase_range = 45:5:135;
    n_phases = length(phase_range);

    energy_curve = zeros(n_phases, 1);
    power_curve = zeros(n_phases, 1);
    efficiency_curve = zeros(n_phases, 1);

    original_phase = params.phaseShift;
    theta = linspace(0, 2*pi, params.simulationPointsPerCycle);

    for i = 1:n_phases
        params.phaseShift = phase_range(i) * pi/180;

        [V_total, V_exp, V_comp, ~, ~] = calc_volumes(theta, params);
        [P, ~, ~] = schmidt_analysis(theta, V_total, V_exp, V_comp, params);
        [W_indicated, P_indicated, ~, ~, ~, efficiency] = ...
            calc_power(P, V_total, theta, params);

        energy_curve(i) = W_indicated;
        power_curve(i) = P_indicated;
        efficiency_curve(i) = efficiency;
    end

    [max_power, max_idx] = max(power_curve);
    optimal_phase = phase_range(max_idx);

    params.phaseShift = original_phase;

    % Return curves with phase angles
    energy_curve = [phase_range(:), energy_curve];
    power_curve = [phase_range(:), power_curve];
    efficiency_curve = [phase_range(:), efficiency_curve];
end