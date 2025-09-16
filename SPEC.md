# Stirling Engine Flywheel Design Project Specification

## Project Overview

### Background
In 1790, Robert Stirling invented a rudimentary yet revolutionary combustion engine known as the Stirling Engine. This project involves analyzing and designing a flywheel for a beta-type Stirling engine using computational modeling techniques.

### Objective
Design and analyze a properly sized flywheel for a beta-type Stirling engine that maintains rotational speed fluctuation within allowable limits throughout the engine cycle.

## Technical Requirements

### Software Requirements
- **MATLAB**: Must be used for all calculations
- **Python**: May be used to wrap MATLAB code for running multiple scripts

### Engineering Assumptions
1. Frictional losses can be neglected
2. All components are massless (except the flywheel)
3. Working fluid behaves as an ideal gas
4. Isothermal expansion and compression
5. Ideal regenerator
6. Lumped, 3-body system
7. Constant pressure throughout the gas
8. Constant temperature in each body

### Flywheel Configuration
- Cylindrical shape with known thickness and width
- Single material with uniform cross-section annulus
- Rotational inertia considered only in outer thickness section
- Spokes have no mass impact on rotational inertia
- Properly sized to maintain speed fluctuation within allowable limits

### Engine Configuration
- Beta-type Stirling engine
- Crank slider mechanism for both displacer and power pistons
- Both cranks connected to same shaft at an angle (Phase Shift)

## Input Parameters

### Power Piston
- Crank Length
- Connecting Rod Length

### Displacer
- Crank Length
- Connecting Rod Length
- Displacer Volume

### Cylinder
- Bore (Diameter)

### Operating Parameters
- Phase Shift
- Compression Ratio
- High Temperature
- Low Temperature
- Gas Pressure at Bottom Dead Center (BDC)
- Atmospheric Pressure
- Regenerator Dead Volume
- Working Fluid

### Flywheel Parameters
- Width (w)
- Outer Diameter (D) - **Unknown, to be determined**
- Rim Thickness (t)
- Coefficient of Fluctuation
- Average Rotational Velocity

## Output Requirements

### Primary Output
- Flywheel diameter calculation

### Required Plots
1. **P-V Diagram**: Pressure vs. specific volume for both Stirling cycle & Stirling engine over complete engine cycle
2. **Torque Analysis**: Torque vs. crank angle for single engine cycle
3. **Speed Variation**: Rotational velocity vs. crank angle for single engine cycle
4. **Optimization**: Energy per cycle vs. phase angle

## Deliverables

### Technical Deliverables
1. **Flywheel Design**: Basic flywheel design with calculated dimensions
2. **Power Analysis**: Power output calculations using two different methods
3. **Visualization**: All four required plots (see Output Requirements)
4. **Analysis Description**: Text description of analysis being performed
5. **Results Summary**: Text summary of results

### Documentation Deliverables

#### MATLAB Script
- Single `.m` file that performs complete analysis
- Generates all required plots automatically when executed

#### Technical Report
**Format Requirements:**
- Length: 6-10 pages
- Style: Executive summary format
- Target Audience: Professors (undergraduate mechanical engineering level)

**Required Content:**
- All four plots with proper labeling and legends
- Key findings (both quantitative and qualitative)
- Analysis methodology
- Results interpretation
- Conclusion addressing:
  - Which design requirements were met
  - Which requirements were missed and explanation

**Structure:**
1. Executive Summary
2. Introduction and Background
3. Methodology
4. Analysis and Results
5. Key Findings
6. Conclusions and Recommendations
7. Appendices (if needed)

## Additional Resources

- **Stirling Engine Overview**: [Wikipedia - Stirling Engine](https://en.wikipedia.org/wiki/Stirling_engine)
- **Beta-Type Configuration**: [Ohio University - Beta Stirling Engines](https://people.ohio.edu/urieli/stirling/engines/beta.html)

## Project Timeline and Submission

*To be determined based on course schedule*

---

*Note: This specification document serves as the primary reference for the Stirling Engine Flywheel Design project. All team members should review and understand these requirements before beginning implementation.*