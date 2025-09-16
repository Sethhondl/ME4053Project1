%% STIRLING ENGINE FLYWHEEL DESIGN ANALYSIS
% Beta-type Stirling Engine with Crank-Slider Mechanism
% This script performs comprehensive analysis including:
%   - Thermodynamic cycle analysis using Schmidt theory
%   - Flywheel sizing based on energy fluctuation
%   - Power output calculation using two methods
%   - Dynamic simulation of speed variation
%   - Phase angle optimization
%   - Generation of all required plots
%
% Author: Stirling Engine Analysis System
% Date: Generated Analysis
% Course: Mechanical Engineering Modeling

clear all; close all; clc;

fprintf('================================================================\n');
fprintf('     STIRLING ENGINE FLYWHEEL DESIGN - ANALYSIS SYSTEM         \n');
fprintf('================================================================\n\n');

%% 1. LOAD PARAMETERS
fprintf('Step 1: Loading Engine Parameters...\n');
fprintf('----------------------------------------\n');
params = engine_parameters();

%% 2. DEFINE SIMULATION PARAMETERS
fprintf('\nStep 2: Setting Up Simulation...\n');
fprintf('----------------------------------------\n');
theta = linspace(0, 2*pi, params.sim.n_points);  % Crank angle array (radians)
theta = theta(:);  % Ensure column vector
fprintf('  Simulation Points: %d per cycle\n', params.sim.n_points);
fprintf('  Angular Resolution: %.2f degrees\n', 360/params.sim.n_points);

%% 3. VOLUME ANALYSIS
fprintf('\nStep 3: Calculating Volumes...\n');
fprintf('----------------------------------------\n');
[V_total, V_exp, V_comp, x_power, x_disp] = calc_volumes(theta, params);

%% 4. THERMODYNAMIC ANALYSIS (Schmidt Theory)
fprintf('\nStep 4: Performing Schmidt Analysis...\n');
fprintf('----------------------------------------\n');
[P, m_total, P_mean] = schmidt_analysis(theta, V_total, V_exp, V_comp, params);

%% 5. WORK AND POWER CALCULATIONS
fprintf('\nStep 5: Calculating Work and Power...\n');
fprintf('----------------------------------------\n');
[W_indicated, P_indicated, W_mep, P_mep, MEP, efficiency] = ...
    calc_power(P, V_total, theta, params);

%% 6. TORQUE ANALYSIS
fprintf('\nStep 6: Analyzing Torque Profile...\n');
fprintf('----------------------------------------\n');
[T_total, T_power, T_disp, T_mean] = calc_torque(P, theta, x_power, x_disp, params);

%% 7. FLYWHEEL DESIGN
fprintf('\nStep 7: Sizing Flywheel...\n');
fprintf('----------------------------------------\n');
[D_outer, I_required, mass_flywheel, energy_fluctuation] = ...
    size_flywheel(T_total, theta, params);

%% 8. DYNAMIC SIMULATION
fprintf('\nStep 8: Simulating Angular Velocity Variation...\n');
fprintf('----------------------------------------\n');
[omega, alpha, rpm, Cs_actual] = simulate_dynamics(T_total, theta, I_required, params);

%% 9. PHASE ANGLE OPTIMIZATION
fprintf('\nStep 9: Optimizing Phase Angle...\n');
fprintf('----------------------------------------\n');
[optimal_phase, energy_curve, power_curve, efficiency_curve] = optimize_phase(params);

%% 10. COMPILE RESULTS STRUCTURE
fprintf('\nStep 10: Compiling Results...\n');
fprintf('----------------------------------------\n');

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

% Power and efficiency
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

% Dynamic simulation
results.omega = omega;
results.alpha = alpha;
results.rpm = rpm;
results.rpm_mean = mean(rpm);
results.Cs_actual = Cs_actual;

% Optimization
results.optimal_phase = optimal_phase;
results.optimization.energy_curve = energy_curve;
results.optimization.power_curve = power_curve;
results.optimization.efficiency_curve = efficiency_curve;

% Angle data
results.theta = theta;

fprintf('  Results structure created successfully\n');

%% 11. GENERATE ALL PLOTS
fprintf('\nStep 11: Generating Plots...\n');
fprintf('----------------------------------------\n');
generate_all_plots(results, params);

%% 12. DISPLAY COMPREHENSIVE RESULTS
fprintf('\nStep 12: Final Results Summary\n');
fprintf('----------------------------------------\n');
display_results(results, params);

%% 13. ANALYSIS COMPLETE
fprintf('================================================================\n');
fprintf('           ANALYSIS COMPLETE - ALL TASKS FINISHED              \n');
fprintf('================================================================\n\n');

fprintf('Summary of Deliverables:\n');
fprintf('  1. Flywheel Design: D = %.3f m, M = %.1f kg\n', ...
        D_outer, mass_flywheel);
fprintf('  2. Power Output: %.2f kW (Method 1), %.2f kW (Method 2)\n', ...
        P_indicated/1000, P_mep/1000);
fprintf('  3. Four Required Plots: Saved to results/ directory\n');
fprintf('  4. Analysis Description: Included in code comments\n');
fprintf('  5. Results Summary: Displayed above and saved to file\n');
fprintf('\nAll project requirements have been addressed.\n');

%% EXPORT KEY VARIABLES TO WORKSPACE
% Make key results available for further analysis
fprintf('\nKey variables exported to workspace:\n');
fprintf('  - params: Engine parameters\n');
fprintf('  - results: Complete results structure\n');
fprintf('  - P, V_total: Pressure and volume arrays\n');
fprintf('  - T_total: Torque array\n');
fprintf('  - rpm: Speed variation array\n');

% End of main script