# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a MATLAB-based Stirling Engine analysis system for designing flywheels for beta-type Stirling engines. The system performs thermodynamic cycle analysis using Schmidt theory, sizes flywheels based on energy fluctuation, and optimizes phase angles for maximum power output.

## Essential Commands

### Running the Analysis
```matlab
% Main analysis - run from Project1 directory
stirling_engine_analysis

% Run with custom parameters (edit engine_parameters.m first)
clear all; stirling_engine_analysis

% Batch mode for debugging (from terminal)
matlab -batch "stirling_engine_analysis"

% Test specific components
params = engine_parameters();
theta = linspace(0, 2*pi, 360);
[V_total, V_exp, V_comp, x_power, x_disp] = calc_volumes(theta, params);
```

### Debugging and Testing
```bash
# Run and capture errors
matlab -batch "try; stirling_engine_analysis; catch ME; fprintf('Error: %s\n', ME.message); disp(ME.stack(1)); end"

# Check key outputs only
matlab -batch "stirling_engine_analysis" | grep -E "Power Output:|Efficiency:|STATUS:"
```

## Architecture Overview

### Calculation Pipeline
The system follows a strict sequential pipeline where each step depends on previous results:

1. **engine_parameters.m** → Defines all configuration (must maintain physical consistency)
2. **calc_volumes.m** → Calculates instantaneous volumes using crank-slider kinematics
3. **schmidt_analysis.m** → Determines pressure from volumes using isothermal assumption
4. **calc_power.m** → Computes work via two methods for validation:
   - Method 1: Numerical integration of P-dV
   - Method 2: Mean Effective Pressure approach
5. **calc_torque.m** → Converts pressure forces to crank torque
6. **size_flywheel.m** → Calculates required inertia from energy fluctuation
7. **simulate_dynamics.m** → Validates speed fluctuation meets Cs requirement
8. **optimize_phase.m** → Tests phase angles 60-120° for optimal power
9. **generate_all_plots.m** → Creates all visualizations with proper ideal cycle
10. **display_results.m** → Validates against requirements and reports

### Critical Implementation Details

#### Beta-Type Engine Volume Model
In `calc_volumes.m`, the beta-type engine has specific characteristics:
- Total volume depends ONLY on power piston position
- Displacer only redistributes gas between hot/cold spaces
- Must use element-wise operations (.*) for array calculations
- Compression space: `V_comp = V_working .* (1 - 0.8 * disp_fraction)`

#### Work Calculation Sign Convention
In `calc_power.m`:
- Clockwise P-V cycle = positive work (power production)
- Counter-clockwise = negative work (power consumption)
- Uses `abs(W_indicated)` to handle sign issues
- Both calculation methods must agree within 5%

#### Schmidt Analysis Assumptions
- Instantaneous pressure equilibrium throughout engine
- Isothermal compression/expansion spaces
- Ideal regenerator with perfect effectiveness
- Mass conservation: `P = (m_total * R) / (V_c/T_c + V_r/T_r + V_h/T_h)`

## Key Parameter Constraints

### Target Specifications
- **Power Output**: 1-10 kW (adjust bore: 120-200mm, pressure: 1-2 MPa)
- **Compression Ratio**: 3-20 (calculated as V_max/V_min)
- **Coefficient of Fluctuation**: ≤ 0.04
- **Efficiency**: Must be < Carnot limit (1 - T_cold/T_hot)
- **Flywheel Diameter**: < 2.0 m practical limit

### Typical Working Values
- Bore: 150-180 mm for 1-10 kW
- Stroke: 60-120 mm
- Operating pressure: 1-1.5 MPa at BDC
- Temperatures: T_hot=600-700K, T_cold=300-350K
- Phase shift: 90° (optimal typically 90-105°)
- RPM: 300-600

## Common Issues and Solutions

### Matrix Operation Errors
```matlab
% WRONG - causes dimension error
V_comp = V_working * fraction;

% CORRECT - element-wise multiplication
V_comp = V_working .* fraction;
```

### Power Outside Target Range
- Too low: Increase bore, stroke, or operating pressure
- Too high: Reduce dimensions or pressure
- Check compression ratio is reasonable (5-15 typical)

### Negative Work/Torque
- Usually indicates P-V cycle going counter-clockwise
- Check volume calculation phase relationships
- Verify displacer leads power piston by phase angle

### Efficiency Issues
- If approaching Carnot: Check heat input calculation
- If too low: Verify temperature differential
- Ensure dead volumes are reasonable (< 20% of swept volume)

## Validation Checklist

After running analysis, verify:
1. STATUS shows "ALL DESIGN REQUIREMENTS MET ✓"
2. Power output: 1-10 kW range
3. Efficiency: 20-40% typical (< Carnot)
4. Compression ratio: 3-20 range
5. Flywheel diameter: < 1.5 m typical
6. Coefficient of fluctuation: ≤ 0.04
7. All four plots generated in results/
8. P-V diagram shows clockwise cycle

## Plot Validation

### P-V Diagram Should Show:
- Clockwise cycle (power producing)
- Ideal cycle: isothermal curves + vertical isochoric lines
- Real cycle: smooth quasi-elliptical shape
- Positive work annotation

### Torque Profile Should Show:
- Oscillating pattern with 2-4 zero crossings
- Peak torques ~10x mean torque
- Mean torque value (may show negative but power is positive)

### Phase Optimization Should Show:
- Maximum power near 90-100°
- Dome-shaped power curve
- Energy and power curves aligned

## Report Generation

For academic reports, the system provides:
- Executive summary format output
- Four publication-ready PNG plots
- Detailed methodology in ANALYSIS_DESCRIPTION.md
- Key findings in RESULTS_SUMMARY.md
- Validation against all requirements