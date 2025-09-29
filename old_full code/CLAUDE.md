# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a MATLAB-based Stirling Engine analysis system for a school project. The PRIMARY GOAL is to **calculate the required flywheel diameter** for a beta-type Stirling engine to maintain speed fluctuation within allowable limits (Cs ≤ coefficient of fluctuation).

**IMPORTANT**: This is an academic exercise focused on flywheel sizing calculations. The goal is NOT to optimize engine performance or achieve specific power output targets, but to correctly size a flywheel based on the given engine parameters.

## Project Specifications (from spec.md)

### Primary Objective
Design and analyze a properly sized flywheel for a beta-type Stirling engine that maintains rotational speed fluctuation within the coefficient of fluctuation limit.

### Key Requirements
1. **Primary Output**: Calculate flywheel outer diameter (D)
2. **Coefficient of Fluctuation**: Must maintain Cs ≤ target value (given in parameters)
3. **Analysis Methods**: Power calculation using TWO different methods for validation
4. **Software**: MATLAB for all calculations

### Required Deliverables
1. **Flywheel Design**: Calculate outer diameter based on energy fluctuation
2. **Power Analysis**: Two calculation methods (P-dV integration and MEP)
3. **Four Required Plots**:
   - P-V Diagram (Stirling cycle vs actual engine)
   - Torque vs crank angle
   - Rotational velocity vs crank angle
   - Energy per cycle vs phase angle (optimization)
4. **Analysis Description**: Text description of methodology
5. **Results Summary**: Executive summary format

## Essential Commands

### Running the Analysis
```matlab
% Main analysis - run from Project1 directory
stirling_engine_analysis

% For clean code version (no comments)
cd clean_code
stirling_engine_analysis
```

## Architecture Overview

### Calculation Pipeline
The system follows a strict sequential pipeline focused on flywheel sizing:

1. **engine_parameters.m** → Defines given engine configuration
2. **calc_volumes.m** → Calculates instantaneous volumes using crank-slider kinematics
3. **schmidt_analysis.m** → Determines pressure using Schmidt theory
4. **calc_power.m** → Computes work via TWO methods (requirement)
5. **calc_torque.m** → Converts pressure forces to crank torque
6. **size_flywheel.m** → **PRIMARY GOAL: Calculates required flywheel diameter**
7. **simulate_dynamics.m** → Validates that Cs requirement is met
8. **optimize_phase.m** → Shows effect of phase angle on performance
9. **generate_all_plots.m** → Creates four required plots
10. **display_results.m** → Presents results in executive summary format

### Critical Implementation Details

#### Flywheel Sizing (PRIMARY FOCUS)
In `size_flywheel.m`:
- Calculate energy fluctuation from torque variation
- Use Cs formula: I_required = ΔE / (Cs * ω²)
- Solve for outer diameter given rim geometry
- Validate that Cs ≤ target value

#### Beta-Type Engine Volume Model
- Total volume depends ONLY on power piston position
- Displacer redistributes gas between hot/cold spaces
- Volume conservation must be maintained

#### Work Calculation (TWO METHODS REQUIRED)
- Method 1: Direct P-dV integration using trapz
- Method 2: Mean Effective Pressure (MEP) approach
- Both must agree within reasonable tolerance (~5%)

## Given Parameters (from engine_parameters.m)

### Current Configuration
- **Power Piston**: Crank 25mm, rod 75mm
- **Displacer**: Crank 20mm, rod 140mm, volume 40 cm³
- **Cylinder**: Bore 50mm
- **Operating**: Phase 90°, CR 1.7, T_hot=900K, T_cold=300K
- **Pressure**: 500 kPa at BDC
- **Flywheel**: Width 25mm, rim 50mm, Cs=0.003, 650 RPM
- **Material**: Steel (7750 kg/m³)

### Engineering Assumptions (from spec.md)
1. Frictional losses neglected
2. Massless components (except flywheel)
3. Ideal gas behavior
4. Isothermal expansion/compression
5. Ideal regenerator
6. Uniform pressure throughout gas
7. Constant temperature in each space

## Common Issues and Solutions

### Focus on Flywheel Sizing
- The power output is NOT a target - it's a calculation result
- The compression ratio of 1.7 is GIVEN, not a requirement to meet
- Focus on calculating correct flywheel diameter for given Cs

### Current Results (Clean Code)
- Power: 0.136 kW (this is fine - not a target)
- Compression Ratio: 1.70 (matches given value)
- Cs: 0.003 (meets requirement ✓)
- Flywheel Diameter: 0.843 m (PRIMARY OUTPUT)
- Efficiency: 2.8% (< Carnot ✓)

### If Parameters Need Adjustment
The given parameters (50mm bore, 500 kPa) produce low power but that's OK for this academic exercise. The focus is on demonstrating proper flywheel sizing methodology, not achieving industrial power levels.

## Validation Checklist

After running analysis, verify:
1. ✓ Flywheel diameter calculated (0.843 m with current params)
2. ✓ Coefficient of fluctuation ≤ target (0.003 ≤ 0.003)
3. ✓ Two power methods agree (100% agreement)
4. ✓ All four plots generated
5. ✓ Efficiency < Carnot limit (2.8% < 66.7%)
6. ✓ Torque profile shows variation
7. ✓ Speed variation within Cs limits

## Plot Requirements

### 1. P-V Diagram
- Ideal Stirling cycle (isothermal + isochoric)
- Actual engine cycle overlay
- Clockwise rotation (power producing)

### 2. Torque Profile
- Torque vs crank angle (0-360°)
- Mean torque line
- Shows torque reversals

### 3. Speed Variation
- RPM vs crank angle
- Mean speed line
- Verifies Cs = (ωmax - ωmin)/ωmean

### 4. Phase Optimization
- Energy/cycle vs phase angle (60-120°)
- Shows optimal phase angle

## Report Focus

The executive summary should emphasize:
1. **Primary Result**: Calculated flywheel diameter = 0.843 m
2. **Validation**: Cs achieved = 0.003 (meets requirement)
3. **Methodology**: Two power methods show 100% agreement
4. **Key Finding**: Flywheel successfully maintains speed stability
5. **Design Success**: All engineering requirements met

## Important Notes

- This is an ACADEMIC PROJECT focused on flywheel sizing methodology
- Power output is a RESULT of the given parameters, not a design target
- The small 50mm bore engine is perfectly valid for demonstrating the analysis
- The goal is correct flywheel sizing, not engine optimization
- All given parameters should be used as-is unless specifically instructed otherwise