# Stirling Engine Flywheel Design Project

## Overview
This MATLAB project implements a comprehensive analysis system for designing a flywheel for a beta-type Stirling engine. The analysis includes thermodynamic cycle calculations, flywheel sizing, power output analysis, and performance optimization.

## Project Structure

```
Project1/
├── stirling_engine_analysis.m    # Main script - RUN THIS
├── engine_parameters.m            # Configuration file (edit parameters here)
├── calc_volumes.m                 # Volume calculations
├── schmidt_analysis.m             # Pressure calculations
├── calc_power.m                   # Power and efficiency calculations
├── calc_torque.m                  # Torque analysis
├── size_flywheel.m               # Flywheel sizing algorithm
├── simulate_dynamics.m            # Speed variation simulation
├── optimize_phase.m              # Phase angle optimization
├── generate_all_plots.m          # Plotting functions
├── display_results.m             # Results display and validation
├── SPEC.md                       # Project specifications
├── ANALYSIS_DESCRIPTION.md       # Detailed methodology
├── RESULTS_SUMMARY.md            # Key findings and conclusions
└── results/                      # Output directory (created automatically)
    ├── pv_diagram.png
    ├── torque_profile.png
    ├── velocity_variation.png
    ├── phase_optimization.png
    └── analysis_summary.txt
```

## How to Run

1. **Open MATLAB** and navigate to the Project1 directory

2. **Run the main script**:
   ```matlab
   stirling_engine_analysis
   ```

3. The analysis will automatically:
   - Load parameters
   - Perform all calculations
   - Generate plots
   - Display results
   - Save outputs to the results/ folder

## Modifying Parameters

To change engine parameters, edit the `engine_parameters.m` file:

### Key Parameters to Adjust:
- **Power Piston**: `params.power.crank_length` (default: 100mm)
- **Displacer**: `params.disp.crank_length` (default: 120mm)
- **Temperatures**: `params.T_hot`, `params.T_cold`
- **Phase Shift**: `params.phase_shift` (default: 90°)
- **Operating Speed**: `params.rpm_avg` (default: 500 RPM)
- **Flywheel Coefficient**: `params.flywheel.coeff_fluctuation` (default: 0.04)

## Output Files

### Plots (PNG format):
1. **pv_diagram.png**: Comparison of ideal Stirling cycle vs actual engine cycle
2. **torque_profile.png**: Torque variation over one complete cycle
3. **velocity_variation.png**: Speed fluctuation with flywheel
4. **phase_optimization.png**: Energy and power vs phase angle

### Text Output:
- **Console Output**: Comprehensive results display with validation
- **analysis_summary.txt**: Key results saved to file

## Deliverables Checklist

✓ **Flywheel Design**: Calculated diameter, mass, and inertia
✓ **Power Output**: Calculated using two independent methods
✓ **Four Required Plots**: All generated and saved
✓ **Analysis Description**: See ANALYSIS_DESCRIPTION.md
✓ **Results Summary**: See RESULTS_SUMMARY.md
✓ **MATLAB Script**: Single script that runs complete analysis

## Typical Results

For the default configuration:
- **Power Output**: 3-5 kW
- **Flywheel Diameter**: 0.5-0.8 m
- **Thermal Efficiency**: 35-40%
- **Speed Fluctuation**: ±2% (Cs = 0.04)

## Validation Features

The system includes automatic validation of:
- Conservation of mass
- Positive pressure throughout cycle
- Efficiency below Carnot limit
- Reasonable flywheel dimensions
- Agreement between power calculation methods

## Troubleshooting

### Common Issues:

1. **"Negative pressure detected"**
   - Check compression ratio is reasonable (2-5)
   - Verify temperature inputs are correct

2. **"Flywheel diameter exceeds limit"**
   - Reduce coefficient of fluctuation requirement
   - Check torque calculations
   - Verify operating speed

3. **"Power output outside target range"**
   - Adjust cylinder dimensions
   - Modify operating temperatures
   - Change phase angle

## Requirements

- MATLAB R2018b or newer
- No additional toolboxes required
- ~50 MB disk space for outputs

## Academic Use

This project is designed for:
- Mechanical Engineering undergraduate courses
- Thermodynamics and modeling education
- Stirling engine design studies

## Report Generation

For the executive summary report:
1. Run the complete analysis
2. Collect the four PNG plots from results/
3. Use RESULTS_SUMMARY.md for key findings
4. Use ANALYSIS_DESCRIPTION.md for methodology
5. Include console output for detailed results

## Contact

For questions about the implementation or methodology, refer to:
- SPEC.md for project requirements
- ANALYSIS_DESCRIPTION.md for technical details
- RESULTS_SUMMARY.md for interpretation guidance

---

*Note: This is an educational tool for Stirling engine analysis. Results should be validated with physical testing for actual engine development.*