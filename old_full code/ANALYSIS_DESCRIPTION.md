# Stirling Engine Flywheel Design - Analysis Methodology

## Executive Summary

This analysis implements a comprehensive computational model for designing a flywheel for a beta-type Stirling engine. The analysis uses MATLAB to perform thermodynamic cycle calculations, flywheel sizing, and performance optimization for a large-scale power generation system (1-10 kW range).

## Analysis Methodology

### 1. Volume Calculations (Kinematic Analysis)

The analysis begins by calculating instantaneous volumes using crank-slider kinematics for both the power piston and displacer:

- **Power Piston Position**: Calculated using crank radius and connecting rod length
- **Displacer Position**: Incorporates phase shift (typically 90°) relative to power piston
- **Volume Spaces**: Tracks compression space (cold), expansion space (hot), and dead volumes

The crank-slider mechanism equations account for the non-linear relationship between crank angle and piston position, providing accurate volume variations throughout the cycle.

### 2. Thermodynamic Analysis (Schmidt Theory)

The pressure calculation employs Schmidt analysis, which assumes:
- Isothermal compression and expansion processes
- Ideal regenerator with perfect heat transfer
- Instantaneous pressure equilibrium throughout the engine
- Ideal gas behavior for the working fluid (air)

The Schmidt equation relates pressure to the distribution of gas mass in different temperature zones:
```
P = (m_total × R) / (V_cold/T_cold + V_regenerator/T_avg + V_hot/T_hot)
```

This provides the pressure variation throughout the cycle, which drives the power output.

### 3. Work and Power Calculation (Dual Methods)

Power output is calculated using two independent methods for validation:

**Method 1 - Direct Integration**:
- Numerically integrates the P-dV diagram using trapezoidal rule
- Provides cycle work from the enclosed area of the P-V diagram
- Power = Work × Frequency

**Method 2 - Mean Effective Pressure (MEP)**:
- Calculates the equivalent constant pressure for the same work output
- MEP = Work / Swept Volume
- Provides verification of Method 1

Both methods should agree within 5% for a properly functioning model.

### 4. Torque Analysis

Torque is calculated from pressure forces acting on the pistons:
- Converts linear piston forces to rotational torque via crank mechanism
- Accounts for connecting rod angle effects
- Combines contributions from both power piston and displacer
- Identifies torque fluctuations that must be smoothed by the flywheel

### 5. Flywheel Sizing Algorithm

The flywheel is sized based on energy fluctuation and coefficient of speed fluctuation:

1. **Energy Fluctuation Calculation**:
   - Integrates torque deviation from mean throughout cycle
   - Identifies maximum energy variation (ΔE)

2. **Required Inertia**:
   - I = ΔE / (Cs × ω²)
   - Where Cs is the coefficient of fluctuation (typically 0.04)

3. **Geometric Sizing**:
   - Solves for outer diameter given material properties
   - Assumes annular (hollow) cylinder configuration
   - Material: Steel (ρ = 7850 kg/m³)

### 6. Dynamic Simulation

Validates the flywheel design by simulating speed variation:
- Solves equation of motion: I×α = T_engine - T_load
- Uses energy-based integration for numerical stability
- Verifies that actual speed fluctuation meets design requirements

### 7. Phase Angle Optimization

Systematically varies the phase angle between pistons (60-120°) to find:
- Optimal angle for maximum power output
- Optimal angle for maximum efficiency
- Trade-offs between power and efficiency

## Key Assumptions and Limitations

### Assumptions:
1. Frictionless operation (no mechanical losses)
2. Perfect sealing (no leakage)
3. Ideal gas behavior
4. Isothermal heat transfer in working spaces
5. Massless components (except flywheel)
6. Constant atmospheric pressure

### Limitations:
1. Schmidt analysis provides approximate pressure (±10-15% typical)
2. Does not account for finite heat transfer rates
3. Neglects flow losses through heat exchangers
4. Simplified displacer torque model

## Validation Approach

The analysis includes multiple validation checks:
1. Conservation of mass throughout cycle
2. Positive work output verification
3. Efficiency below Carnot limit
4. Physical reasonableness of dimensions
5. Agreement between two power calculation methods

## Performance Metrics

The analysis provides comprehensive performance metrics:
- Power output (kW)
- Thermal efficiency (% and % of Carnot)
- Flywheel specifications (diameter, mass, inertia)
- Speed fluctuation characteristics
- Optimal operating parameters

## Computational Efficiency

The MATLAB implementation uses:
- Vectorized operations for speed
- Modular functions for maintainability
- Clear variable naming for readability
- Comprehensive error checking
- Professional plotting with proper labels

## Applicability

This analysis is suitable for:
- Initial design of Stirling engine systems
- Flywheel sizing for various operating conditions
- Performance optimization studies
- Educational demonstrations of thermodynamic principles
- Comparison of different design parameters

The modular structure allows easy modification of parameters to explore different engine configurations and operating conditions.