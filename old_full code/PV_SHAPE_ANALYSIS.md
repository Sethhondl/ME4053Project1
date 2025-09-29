# P-V Diagram Shape Analysis: Airfoil vs Kidney Bean

## Current Situation

The P-V diagram shows an elongated "airfoil" shape rather than the typical rounded "kidney bean" shape expected for Stirling engines.

## Root Cause: Compression Ratio Limited to 1.7

The compression ratio of 1.7 (specified in givenPar.csv) is the primary cause of the airfoil shape. This low CR results from:

### Given Constraints:
- Power piston stroke: 50mm (swept volume = 98.2 mL)
- Required CR: 1.7
- This forces total dead volume = swept volume / 0.7 = 140.2 mL

### Dead Volume Breakdown:
- Hot space dead: 48 mL
- Cold space dead: 72 mL
- Regenerator: 20 mL
- **Total: 140 mL** (143% of swept volume!)

## Why This Creates an Airfoil Shape

1. **Limited Volume Range**: With CR=1.7, volume only varies from 140 to 238 mL (70% variation)

2. **Limited Pressure Range**: Pressure only varies by factor of 1.88 (450 to 850 kPa)

3. **Large Dead Volume Buffer**: The 140 mL dead volume acts as a "buffer" that:
   - Dampens pressure swings
   - Reduces the rate of pressure change
   - Creates gradual transitions rather than sharp corners

## Comparison with Different Dead Volumes

| Configuration | CR | Dead/Swept | Work/Cycle | Shape |
|--------------|-----|------------|------------|--------|
| Current (Given) | 1.7 | 143% | 13.0 J | Airfoil |
| Reduced Dead | 3.1 | 47% | 33.1 J | Rounded |
| Minimal Dead | 7.1 | 16% | 66.2 J | Kidney Bean |

## Typical Stirling Engine Values

Real Stirling engines typically have:
- Compression ratio: 2.5 to 4.0
- Dead volume: 20-50% of swept volume
- Result: Kidney bean shaped P-V diagram

## Conclusion

The airfoil shape is **correct** for the given parameters (CR=1.7), but these parameters are not typical of optimized Stirling engines. The constraint of CR=1.7 forces a large dead volume that fundamentally limits the P-V diagram shape.

To achieve a kidney bean shape would require:
1. Increasing the compression ratio to 2.5-4.0
2. Reducing dead volumes to 20-50% of swept volume
3. But this would violate the given CR=1.7 requirement

## Impact on Performance

The airfoil shape (low CR) significantly reduces:
- Work per cycle: Only 13 J vs potential 66 J
- Power output: Limited to 0.14 kW
- Efficiency: Reduced thermodynamic efficiency

The given parameters appear to represent a heavily constrained or de-rated engine design.