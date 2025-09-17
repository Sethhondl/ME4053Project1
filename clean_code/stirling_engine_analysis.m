% Stirling Engine Flywheel Design Analysis
clear all; close all; clc;

% Load parameters
params = engine_parameters();

% Setup simulation
theta = linspace(0, 2*pi, params.simulationPointsPerCycle);
theta = theta(:);

% Volume analysis
[V_total, V_exp, V_comp, x_power, x_disp] = calc_volumes(theta, params);

% Pressure analysis (Schmidt theory)
[P, m_total, P_mean] = schmidt_analysis(theta, V_total, V_exp, V_comp, params);

% Power calculations
[W_indicated, P_indicated, W_mep, P_mep, MEP, efficiency] = ...
    calc_power(P, V_total, theta, params);

% Torque analysis
[T_total, T_power, T_disp, T_mean] = calc_torque(P, theta, x_power, x_disp, params);

% Flywheel design
[D_outer, I_required, mass_flywheel, energy_fluctuation] = ...
    size_flywheel(T_total, theta, params);

% Dynamic simulation
[omega, alpha, rpm, Cs_actual] = simulate_dynamics(T_total, theta, I_required, params);

% Phase optimization
[optimal_phase, energy_curve, power_curve, efficiency_curve] = optimize_phase(params);

% Compile results
results = struct();
results.V_total = V_total;
results.V_exp = V_exp;
results.V_comp = V_comp;
results.x_power = x_power;
results.x_disp = x_disp;
results.P = P;
results.m_total = m_total;
results.P_mean = P_mean;
results.W_indicated = W_indicated;
results.P_indicated = P_indicated;
results.W_mep = W_mep;
results.P_mep = P_mep;
results.MEP = MEP;
results.efficiency = efficiency;
results.T_total = T_total;
results.T_power = T_power;
results.T_disp = T_disp;
results.T_mean = T_mean;
results.flywheel.D_outer = D_outer;
results.flywheel.I_required = I_required;
results.flywheel.mass = mass_flywheel;
results.flywheel.energy_fluctuation = energy_fluctuation;
results.omega = omega;
results.alpha = alpha;
results.rpm = rpm;
results.rpm_mean = mean(rpm);
results.Cs_actual = Cs_actual;
results.optimal_phase = optimal_phase;
results.optimization.energy_curve = energy_curve;
results.optimization.power_curve = power_curve;
results.optimization.efficiency_curve = efficiency_curve;
results.theta = theta;

% Generate plots
generate_all_plots(results, params);

% Display results
display_results(results, params);