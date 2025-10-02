# Appendix: Artificial Intelligence Usage in Project Development

## Overview

This appendix documents the use of artificial intelligence (AI) tools throughout the development of this Stirling engine flywheel sizing analysis project. In adherence to academic integrity standards, we provide full transparency regarding how AI assisted in code development, analysis, visualization, and documentation.

## AI Tools Utilized

### Team Member 1: Claude Code by Anthropic
- **Tool**: Claude Code (claude.ai/code)
- **Model**: Claude Sonnet 4.5
- **Primary Responsibilities**:
  - MATLAB code development and implementation
  - Mathematical modeling and Schmidt analysis
  - Data visualization and plot generation
  - Git version control and documentation

### Team Member 2: Cursor AI with GPT-5
- **Tools**: Cursor IDE with OpenAI GPT-5
- **Primary Responsibilities**:
  - Code review and validation
  - Additional code sections and modules
  - LaTeX report formatting and typesetting
  - Documentation refinement

## Detailed Applications

### 1. MATLAB Code Development

AI assistance was instrumental in developing the computational framework for this project:

#### Core Calculations
- **Flywheel Sizing Algorithm**: AI helped implement the energy fluctuation method and coefficient of fluctuation calculations (Cs = ΔE / (Cs × ω²))
- **Schmidt Analysis**: Assisted in translating thermodynamic theory into MATLAB code for pressure and temperature calculations
- **Volume Calculations**: Implemented crank-slider kinematics for beta-type Stirling engine geometry
- **Torque Analysis**: Developed algorithms to convert pressure forces into instantaneous crank torque

#### Code Quality and Structure
- Refactored variable names to descriptive camelCase for readability
- Implemented modular function architecture for maintainability
- Created comprehensive documentation within code
- Validated two independent power calculation methods (P-dV integration and MEP)

### 2. Data Visualization and Plot Enhancement

Extensive AI collaboration produced professional-quality engineering visualizations:

#### Plot Improvements
- **Annotation System**: Added key measurements and labels to all plots:
  - Stroke measurements on piston position plots
  - V_max, V_min, and compression ratio labels
  - P_max, P_min, and pressure ratio annotations
  - Torque extremes and zero crossings
  - RPM extremes and coefficient of fluctuation display

- **P-v Diagram Redesign**: Created comprehensive pressure-volume diagram with:
  - Ideal Stirling cycle overlay (isothermal and isochoric processes)
  - State labels (1-4) at cycle corners
  - Process names with temperature indicators
  - Proper cycle direction indication

- **Label Positioning**: Iterative refinement to ensure:
  - No overlap between labels and plot lines
  - Optimal use of white space
  - Consistent formatting across all figures
  - Professional appearance suitable for engineering reports

- **Technical Specifications**:
  - High-resolution PNG exports (300-600 DPI)
  - Consistent axis ranges (0-360° for crank angle plots)
  - Enhanced line widths (2.5pt) for clarity
  - Optimized legend placement using MATLAB's 'best' algorithm

### 3. Code Review and Validation

AI tools assisted in ensuring code correctness:

- **Algorithm Verification**: Cross-checked thermodynamic calculations against theory
- **Numerical Validation**: Verified that:
  - Two power methods agreed to 100%
  - Efficiency remained below Carnot limit (2.8% < 66.7%)
  - Coefficient of fluctuation met design requirement (0.003 ≤ 0.003)
  - Volume conservation maintained throughout cycle

- **Edge Case Testing**: Identified and resolved numerical stability issues
- **Documentation Review**: Ensured code comments accurately reflected implementation

### 4. Report Formatting and Documentation

AI assisted in producing professional documentation:

- **LaTeX Formatting**: Structured technical report with proper:
  - Equation formatting and numbering
  - Figure and table placement
  - Citation management
  - Section organization

- **Project Documentation**: Created comprehensive CLAUDE.md file documenting:
  - Project architecture and pipeline
  - Critical implementation details
  - Validation checklist
  - Troubleshooting guidance

- **Git Integration**: Maintained clear commit history with descriptive messages

## Rationale for AI Usage

### Engineering Productivity
1. **Rapid Prototyping**: AI enabled quick iteration on algorithms and visualizations
2. **Code Quality**: AI suggestions improved code structure and readability
3. **Documentation**: AI helped maintain comprehensive inline documentation
4. **Debugging Efficiency**: AI accelerated identification and resolution of errors

### Educational Value
1. **Learning Reinforcement**: AI explanations deepened understanding of concepts
2. **Best Practices**: Exposure to professional coding standards and conventions
3. **Iterative Refinement**: AI feedback encouraged multiple improvement cycles
4. **Time Management**: More time available for conceptual understanding vs. syntax debugging

### Professional Standards
1. **Version Control**: AI facilitated proper Git workflow with meaningful commits
2. **Documentation**: Produced industry-standard code documentation
3. **Visualization**: Created publication-quality engineering figures
4. **Report Quality**: Professional LaTeX formatting suitable for technical audiences

## Academic Integrity and Methodology

### Student Oversight and Understanding
- **Active Direction**: Students directed all AI interactions with specific goals and requirements
- **Critical Review**: All AI-generated code was reviewed, tested, and validated by students
- **Conceptual Understanding**: AI served as a tool to implement student understanding, not replace it
- **Problem Solving**: Students made all engineering decisions regarding approach and methodology

### Transparency Measures
1. **Git History**: Complete commit history available showing incremental development
2. **AI Attribution**: Git commits include "Generated with Claude Code" attribution
3. **Documentation**: This appendix provides full disclosure of AI usage
4. **Code Comments**: Inline documentation shows reasoning and approach

### Verification of Results
All computational results were validated through:
- Comparison with theoretical expectations (e.g., Carnot efficiency limit)
- Two independent power calculation methods showing 100% agreement
- Physical reasonableness checks on all outputs
- Manual spot-checks of critical calculations

## Specific Examples of AI Contributions

### Example 1: Plot Annotation Refinement
**Task**: Position labels on plots without overlapping plot lines

**AI Process**:
1. Initial attempt placed labels using automatic positioning
2. User feedback identified overlaps
3. AI iteratively adjusted label coordinates based on visual inspection
4. Final positions achieved through 5-7 refinement cycles

**Result**: Clean, professional plots with all information clearly visible

### Example 2: Flywheel Sizing Implementation
**Task**: Calculate required flywheel diameter for target coefficient of fluctuation

**AI Process**:
1. Student provided theoretical background and equations
2. AI translated equations into MATLAB code
3. Student reviewed implementation for correctness
4. AI refined code based on numerical validation
5. Final implementation verified against hand calculations

**Result**: Correct flywheel diameter (0.843 m) meeting all design requirements

### Example 3: Schmidt Analysis Validation
**Task**: Ensure thermodynamic calculations followed Schmidt theory correctly

**AI Process**:
1. AI implemented initial pressure calculations
2. Student identified theoretical inconsistencies
3. AI revised calculations based on student corrections
4. Multiple validation checks confirmed accuracy
5. Student verified results against literature values

**Result**: Thermodynamically consistent pressure and temperature profiles

## Conclusion

AI tools served as powerful assistants in this project, enabling higher quality code, better visualizations, and more comprehensive documentation than would have been feasible within the project timeline using traditional methods alone. Critically, AI did not replace student engineering judgment or understanding—rather, it amplified student capability by handling repetitive tasks, suggesting improvements, and facilitating rapid iteration.

The transparency provided in this appendix, combined with the complete Git history and comprehensive documentation, demonstrates responsible AI usage that enhances rather than diminishes the educational value of this project. All engineering decisions, conceptual approaches, and validation steps were directed by students who maintain full understanding and ownership of the work.

## References

- Anthropic. (2024). Claude Code. https://claude.ai/code
- OpenAI. (2024). GPT-5 via Cursor IDE
- Git repository: https://github.com/Sethhondl/ME4053Project1

---

*This appendix was prepared with AI assistance (Claude Code) to ensure comprehensive and accurate documentation of AI usage throughout the project, in accordance with academic integrity policies.*
