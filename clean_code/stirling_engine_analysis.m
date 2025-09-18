%% STIRLING ENGINE FLYWHEEL DESIGN ANALYSIS
% Complete thermodynamic and dynamic analysis for beta-type Stirling engine
% Includes flywheel sizing, phase optimization, and performance validation

clear all; close all; clc;

%% =============== INITIALIZATION ===============

% Load engine configuration
params = engine_parameters();

% Create angular position array
theta = linspace(0, 2*pi, params.simulationPointsPerCycle)';


%% =============== THERMODYNAMIC ANALYSIS ===============

% Step 1: Calculate instantaneous volumes
[V_total, V_exp, V_comp, x_power, x_disp] = calc_volumes(theta, params);

% Step 2: Apply Schmidt analysis for pressure
[P, m_total, P_mean] = schmidt_analysis(theta, V_total, V_exp, V_comp, params);

% Step 3: Calculate work and power output
[W_indicated, P_indicated, W_mep, P_mep, MEP, efficiency] = ...
    calc_power(P, V_total, theta, params);


%% =============== MECHANICAL ANALYSIS ===============

% Step 4: Calculate torque profile
[T_total, T_power, T_disp, T_mean] = calc_torque(P, theta, x_power, x_disp, params);

% Step 5: Size flywheel based on energy fluctuation
[D_outer, I_required, mass_flywheel, energy_fluctuation] = ...
    size_flywheel(T_total, theta, params);

% Step 6: Simulate dynamic speed variation
[omega, alpha, rpm, Cs_actual] = simulate_dynamics(T_total, theta, I_required, params);


%% =============== OPTIMIZATION ===============

% Step 7: Find optimal phase angle
[optimal_phase, energy_curve, power_curve, efficiency_curve] = optimize_phase(params);


%% =============== COMPILE RESULTS ===============

results = struct();

% Volume data
results.V_total = V_total;
results.V_exp = V_exp;
results.V_comp = V_comp;
results.x_power = x_power;
results.x_disp = x_disp;

% Pressure data
results.P = P;
results.m_total = m_total;
results.P_mean = P_mean;

% Power data
results.W_indicated = W_indicated;
results.P_indicated = P_indicated;
results.W_mep = W_mep;
results.P_mep = P_mep;
results.MEP = MEP;
results.efficiency = efficiency;

% Torque data
results.T_total = T_total;
results.T_power = T_power;
results.T_disp = T_disp;
results.T_mean = T_mean;

% Flywheel design
results.flywheel.D_outer = D_outer;
results.flywheel.I_required = I_required;
results.flywheel.mass = mass_flywheel;
results.flywheel.energy_fluctuation = energy_fluctuation;

% Dynamic performance
results.omega = omega;
results.alpha = alpha;
results.rpm = rpm;
results.rpm_mean = mean(rpm);
results.Cs_actual = Cs_actual;

% Optimization results
results.optimal_phase = optimal_phase;
results.optimization.energy_curve = energy_curve;
results.optimization.power_curve = power_curve;
results.optimization.efficiency_curve = efficiency_curve;

% Simulation data
results.theta = theta;


%% =============== OUTPUT GENERATION ===============

% Generate all visualization plots
generate_all_plots(results, params);

% Display comprehensive results summary
display_results(results, params);