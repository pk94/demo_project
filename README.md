# Roll Damping Derivative (C_l_p) - Theory and Formulas

## Overview

The **roll damping derivative** $C_{l_p}$ quantifies how the rolling moment coefficient changes with roll angular velocity. This document covers **subsonic, transonic, and supersonic** flight regimes.

---

## Quick Start

### Python
```bash
python roll_damping_derivative.py
```

### Fortran
```bash
# Using gfortran
gfortran -O2 -o roll_damping roll_damping_derivative.f90
./roll_damping

# Using Intel Fortran
ifort -O2 -o roll_damping roll_damping_derivative.f90
./roll_damping

# Windows (gfortran)
gfortran -O2 -o roll_damping.exe roll_damping_derivative.f90
roll_damping.exe
```

---

## 1. Definition

$$
C_{l_p} = \frac{\partial C_l}{\partial \hat{p}}
$$

Where the **non-dimensional roll rate** is:

$$
\hat{p} = \frac{p \cdot b}{2V}
$$

| Symbol | Description | Units |
|--------|-------------|-------|
| $C_l$ | Rolling moment coefficient | dimensionless |
| $p$ | Roll angular velocity | rad/s |
| $b$ | Wingspan | m |
| $V$ | True airspeed | m/s |

---

## 2. Flow Regime Classification

| Regime | Mach Range | Characteristics |
|--------|------------|-----------------|
| **Subsonic** | $M \lt 0.8$ | Prandtl-Glauert compressibility |
| **Transonic** | $0.8 \leq M \leq 1.2$ | Mixed sub/supersonic flow, shocks |
| **Supersonic** | $1.2 \lt M \lt 5.0$ | Ackeret linear theory applies |
| **Hypersonic** | $M \geq 5.0$ | Newtonian flow, high temperature |

---

## 3. Compressibility Factor (β)

The compressibility factor differs between regimes:

**Subsonic (Prandtl-Glauert):**
$$
\boxed{\beta = \sqrt{1 - M^2}}
$$

**Supersonic (Ackeret):**
$$
\boxed{\beta = \sqrt{M^2 - 1}}
$$

---

## 4. Mach Cone and Leading Edge Classification

### 4.1 Mach Cone Angle

In supersonic flow, disturbances propagate within a cone:

$$
\boxed{\mu = \arcsin\left(\frac{1}{M}\right)}
$$

```
                    Mach Cone
                      /
        Freestream   /  μ = Mach angle
        ─────────────●────────
        M > 1        \
                      \
```

| Mach | μ (degrees) |
|------|-------------|
| 1.0 | 90° |
| 1.2 | 56.4° |
| 1.5 | 41.8° |
| 2.0 | 30.0° |
| 3.0 | 19.5° |

### 4.2 Leading Edge Classification

For supersonic flight, wing leading edge behavior depends on sweep:

$$
\text{Subsonic LE: } \Lambda_{LE} > \mu
$$
$$
\text{Supersonic LE: } \Lambda_{LE} < \mu
$$

```
      Subsonic LE                    Supersonic LE
    (High sweep > μ)               (Low sweep < μ)
         ╱                              │
        ╱  Λ > μ                        │  Λ < μ
       ╱                                │
    Flow "sees" subsonic           Flow "sees" supersonic
    component normal to LE         normal to LE
```

**Why it matters:**
- **Subsonic LE**: Better roll damping, more lift at given Mach
- **Supersonic LE**: Reduced aerodynamic effectiveness

---

## 5. Lift Curve Slope Formulas

### 5.1 Subsonic - 2D Airfoil

**Thin Airfoil Theory with Compressibility:**
$$
\boxed{C_{l_\alpha,2D} = \frac{2\pi}{\beta} = \frac{2\pi}{\sqrt{1-M^2}}}
$$

### 5.2 Supersonic - 2D Airfoil (Ackeret Theory)

$$
\boxed{C_{l_\alpha,2D} = \frac{4}{\beta} = \frac{4}{\sqrt{M^2-1}}}
$$

**Key difference**: Supersonic lift curve slope is **lower** than subsonic.

| Mach | $C_{l_\alpha}$ Subsonic | $C_{l_\alpha}$ Supersonic |
|------|-------------------------|---------------------------|
| 0.5 | 7.26 | - |
| 0.8 | 10.47 | - |
| 1.2 | - | 6.03 |
| 2.0 | - | 2.31 |
| 3.0 | - | 1.41 |

### 5.3 Subsonic - 3D Wing (Helmbold Equation)

$$
\boxed{C_{L_\alpha,3D} = \frac{C_{l_\alpha,2D} \cdot AR_{eff}}{2 + \sqrt{4 + \left(\frac{AR_{eff} \cdot \beta \cdot \cos\Lambda}{C_{l_\alpha,2D}}\right)^2}}}
$$

Where:
- $AR_{eff} = AR \cdot \cos\Lambda$

### 5.4 Supersonic - 3D Wing

**For Supersonic Leading Edge:**
$$
\boxed{C_{L_\alpha,3D} = \frac{4}{\beta}\left(1 - \frac{1}{AR \cdot \beta}\right)}
$$

**For Subsonic Leading Edge (swept wing):**

The flow normal to the leading edge is analyzed:

$$
M_n = M \cdot \cos\Lambda_{LE}
$$

If $M_n \lt 1$:
$$
C_{L_\alpha,3D} = \frac{2\pi \cdot AR}{2 + \sqrt{AR^2 \cdot \beta_n^2 + 4}}
$$

where $\beta_n = \sqrt{1 - M_n^2}$

---

## 6. Roll Damping Derivative Formulas

### 6.1 Subsonic Flow

**Tapered Wing Formula:**
$$
\boxed{C_{l_p} = -\frac{C_{L_\alpha}}{12} \cdot \frac{1 + 3\lambda}{1 + \lambda} \cdot \cos^2\Lambda}
$$

**Simple Rectangular Wing:**
$$
C_{l_p} \approx -\frac{C_{L_\alpha}}{4}
$$

### 6.2 Supersonic Flow - Analytical Method

**For Supersonic Leading Edge:**
$$
\boxed{C_{l_p} = -\frac{4}{3\beta} \cdot \frac{AR}{AR + \frac{4}{\beta}} \cdot \frac{1 + 2\lambda}{1 + \lambda}}
$$

**For Subsonic Leading Edge:**

When $M_n = M \cos\Lambda_{LE} \lt 1$ (normal Mach subsonic):
$$
C_{l_p} = -\frac{2\pi}{3(AR + 2/\beta_n)} \cdot \frac{1 + 2\lambda}{1 + \lambda} \cdot \cos^2\Lambda
$$

When $M_n \gt 1$ (normal Mach supersonic):
$$
C_{l_p} = -\frac{4}{3\beta_n} \cdot \frac{AR}{AR + 4/\beta_n} \cdot \frac{1 + 2\lambda}{1 + \lambda} \cdot \cos^2\Lambda
$$

### 6.3 Supersonic Flow - Strip Theory

Numerical integration along span:
$$
\boxed{C_{l_p} = -\frac{4}{S \cdot b^2} \int_0^{b/2} c(y) \cdot C_{L_\alpha}(y) \cdot y^2 \, dy}
$$

With **supersonic local lift curve slope**:
$$
C_{L_\alpha}(y) = \frac{4}{\beta} \cdot \left(1 - \left(\frac{2y}{b}\right)^2\right)
$$

The $(1 - (2y/b)^2)$ term is a tip loss factor for finite span.

### 6.4 Transonic Flow (Interpolation)

Since transonic flow is complex with shock formation:
$$
C_{l_p,transonic} = C_{l_p,sub}(1-t) + C_{l_p,sup} \cdot t
$$

Where:
$$
t = \frac{M - 0.8}{0.4}, \quad t_{smooth} = 3t^2 - 2t^3
$$

### 6.5 Hypersonic Flow

For $M \gt 5$, use **Mach number independence** with Newtonian correction:
$$
\boxed{C_{l_p,hyper} = C_{l_p}(M=5) \cdot \sqrt{\frac{5}{M}}}
$$

---

## 7. Physical Interpretation

### Roll Damping Mechanism

```
                    Roll Rate p (clockwise viewed from behind)
                              ↻

      Upgoing Wing                      Downgoing Wing
      ┌──────────────┐                  ┌──────────────┐
      │  Δα < 0      │                  │  Δα > 0      │
      │  ΔL < 0      │                  │  ΔL > 0      │
      └──────────────┘                  └──────────────┘
            ↑                                 ↓
        Less Lift                         More Lift

        Creates restoring moment opposing roll → DAMPING
```

**Local angle of attack change:**
$$
\Delta\alpha(y) = \frac{p \cdot y}{V} = \frac{2y}{b} \cdot \hat{p}
$$

---

## 8. Effect of Mach Number on C_l_p

```
|C_l_p|
   ^
   │
0.5│    ●────●
   │          ╲
0.4│           ╲  Transonic
   │            ╲   dip
0.3│             ●
   │            ╱ ╲
0.2│           ╱   ●────●────●
   │          ╱        Supersonic
0.1│         ╱
   │
   └──────────────────────────────> Mach
       0.5   0.8  1.0  1.5   2.0  2.5
             │    │
         Subsonic Supersonic
```

**Key observations:**
1. Roll damping **decreases** (|C_l_p| reduces) with increasing Mach
2. **Transonic dip**: Minimum damping around M ≈ 1.0
3. Supersonic: Damping stabilizes at lower values
4. Swept wings provide better damping in supersonic regime

---

## 9. Sweep Angle Effects

| Sweep | Subsonic Effect | Supersonic Effect |
|-------|-----------------|-------------------|
| $\Lambda = 0°$ | Max damping | Supersonic LE, reduced damping |
| $\Lambda = 30°$ | Moderate reduction | Transition region |
| $\Lambda = 45°$ | $\cos^2$ reduction | Subsonic LE at M < 1.4, better damping |
| $\Lambda = 60°$ | Significant reduction | Subsonic LE at M < 2.0, best supersonic damping |

**Sweep angle conversion (LE to quarter-chord):**
$$
\tan\Lambda_{c/4} = \tan\Lambda_{LE} - \frac{1-\lambda}{AR(1+\lambda)}
$$

---

## 10. Typical Values

### Subsonic Aircraft
| Aircraft Type | Mach | Typical $C_{l_p}$ |
|---------------|------|-------------------|
| General Aviation | 0.2 | -0.45 to -0.50 |
| Transport | 0.6 | -0.35 to -0.45 |
| Business Jet | 0.8 | -0.30 to -0.40 |

### Supersonic Aircraft
| Aircraft Type | Mach | Typical $C_{l_p}$ |
|---------------|------|-------------------|
| Fighter (swept) | 1.2 | -0.20 to -0.30 |
| Fighter (delta) | 1.6 | -0.15 to -0.25 |
| Fighter (high sweep) | 2.0 | -0.10 to -0.20 |

---

## 11. Summary Flowchart

```
                    ┌─────────────────┐
                    │   Input: M, Wing│
                    └────────┬────────┘
                             │
                             ▼
                    ┌─────────────────┐
                    │ Classify Regime │
                    │ Sub/Trans/Super │
                    └────────┬────────┘
                             │
         ┌───────────────────┼───────────────────┐
         │                   │                   │
         ▼                   ▼                   ▼
   ┌──────────┐       ┌──────────┐       ┌──────────────┐
   │ SUBSONIC │       │TRANSONIC │       │ SUPERSONIC   │
   │ M < 0.8  │       │0.8≤M≤1.2 │       │ M > 1.2      │
   └────┬─────┘       └────┬─────┘       └──────┬───────┘
        │                  │                    │
        ▼                  ▼                    ▼
   ┌──────────┐       ┌──────────┐       ┌──────────────┐
   │β=√(1-M²) │       │Interpolate│      │β = √(M²-1)   │
   │Helmbold  │       │ Sub↔Super │      │Check LE type │
   │C_Lα/12   │       │           │      │Ackeret theory│
   └────┬─────┘       └────┬─────┘       └──────┬───────┘
        │                  │                    │
        └──────────────────┴────────────────────┘
                           │
                           ▼
                    ┌─────────────────┐
                    │      C_l_p      │
                    │    (1/rad)      │
                    └─────────────────┘
```

---

## 12. Dimensional Conversion

To convert to dimensional roll damping derivative:

$$
\boxed{L_p = \frac{\partial L}{\partial p} = C_{l_p} \cdot \frac{\rho V S b^2}{4}}
$$

| Symbol | Description | Units |
|--------|-------------|-------|
| $L_p$ | Dimensional roll damping | N·m·s/rad |
| $\rho$ | Air density | kg/m³ |
| $V$ | True airspeed | m/s |
| $S$ | Wing area | m² |
| $b$ | Wingspan | m |

---

## 13. Standard Atmosphere (ISA)

For altitude-dependent calculations:

**Troposphere (h < 11 km):**
$$
T = 288.15 - 0.0065h \quad \text{(K)}
$$
$$
\rho = \rho_0 \left(\frac{T}{288.15}\right)^{4.2561}
$$

**Stratosphere (11 km < h < 20 km):**
$$
T = 216.65 \text{ K (constant)}
$$
$$
\rho = \rho_{11km} \cdot \exp\left(-\frac{g(h-11000)}{RT}\right)
$$

---

## 14. Code Structure

### Python (`roll_damping_derivative.py`)

```
Classes:
├── FlowRegime (Enum)      - Flow classification
└── WingGeometry           - Wing parameters dataclass

Functions:
├── Compressibility
│   ├── compressibility_factor(mach)
│   ├── mach_cone_angle(mach)
│   └── is_wing_subsonic_le(mach, sweep)
├── Lift Curve Slope
│   ├── lift_curve_slope_2d_subsonic(mach)
│   ├── lift_curve_slope_2d_supersonic(mach)
│   ├── lift_curve_slope_3d_subsonic(AR, sweep, mach)
│   ├── lift_curve_slope_3d_supersonic(wing, mach)
│   └── lift_curve_slope_3d(wing, mach)          # Unified
├── Roll Damping
│   ├── calculate_clp_subsonic(wing, mach)
│   ├── calculate_clp_supersonic_analytical(wing, mach)
│   ├── calculate_clp_supersonic_strip(wing, mach)
│   └── calculate_clp(wing, mach)                # Unified
└── Utilities
    ├── dimensional_roll_damping(clp, rho, V, wing)
    └── standard_atmosphere(altitude)
```

### Fortran (`roll_damping_derivative.f90`)

```
Modules:
├── constants              - PI, DEG2RAD, regime parameters
├── wing_geometry_mod      - WingGeometry type, geometry functions
├── flow_regime_mod        - Regime classification, Mach cone
├── lift_curve_slope_mod   - All C_L_alpha calculations
├── roll_damping_mod       - All C_l_p calculations
└── atmosphere_mod         - ISA standard atmosphere

Main Program:
└── roll_damping_calculator - Example calculations
```

### Key Differences Python vs Fortran

| Feature | Python | Fortran |
|---------|--------|---------|
| Precision | float64 (default) | `real(dp)` (selected_real_kind) |
| Classes | `@dataclass` | `type :: WingGeometry` |
| Enums | `Enum` class | `integer, parameter` |
| Array ops | NumPy | Do loops |
| Integration | `np.trapezoid` | Manual trapezoidal |

---

## 15. Usage Examples

### Python - Custom Wing Analysis

```python
from roll_damping_derivative import WingGeometry, calculate_clp
import numpy as np

# Define custom wing
my_wing = WingGeometry(
    span=12.0,
    root_chord=3.0,
    tip_chord=1.0,
    sweep_leading_edge=np.radians(35)
)

# Calculate at Mach 1.8
clp, regime = calculate_clp(my_wing, mach=1.8)
print(f"C_l_p = {clp:.4f} ({regime.value})")
```

### Fortran - Custom Wing Analysis

```fortran
program custom_analysis
    use constants
    use wing_geometry_mod
    use roll_damping_mod
    implicit none

    type(WingGeometry) :: my_wing
    real(dp) :: clp
    integer :: regime

    my_wing = WingGeometry( &
        span = 12.0_dp, &
        root_chord = 3.0_dp, &
        tip_chord = 1.0_dp, &
        sweep_le = 35.0_dp * DEG2RAD &
    )

    call calculate_clp(my_wing, 1.8_dp, clp, regime)
    write(*,'(A,F8.4)') "C_l_p = ", clp

end program custom_analysis
```

---

## References

1. **Ackeret, J.** - "Luftkräfte auf Flügel" (1925) - Supersonic airfoil theory
2. **Prandtl, L. & Glauert, H.** - Compressibility corrections
3. **Jones, R.T.** - "Properties of Low-Aspect-Ratio Wings at Speeds Below and Above the Speed of Sound" NACA TR-835
4. **USAF DATCOM** - Data Compendium for stability derivatives
5. **Nelson, R.C.** - *Flight Stability and Automatic Control*
6. **Etkin, B. & Reid, L.D.** - *Dynamics of Flight: Stability and Control*
7. **Ashley, H. & Landahl, M.** - *Aerodynamics of Wings and Bodies* (supersonic theory)
