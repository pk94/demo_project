"""
Roll Damping Derivative Calculator (C_l_p)
==========================================
Calculates the derivative of rolling moment coefficient with respect to
non-dimensional roll angular velocity for aircraft wings.

Supports: SUBSONIC, TRANSONIC, and SUPERSONIC flight regimes.

In flight dynamics notation:
    C_l_p = ∂C_l / ∂p̂

where p̂ = pb/(2V) is the non-dimensional roll rate
    p = roll angular velocity (rad/s)
    b = wingspan (m)
    V = true airspeed (m/s)
"""

import numpy as np
from dataclasses import dataclass
from typing import Tuple
from enum import Enum


class FlowRegime(Enum):
    """Flow regime classification based on Mach number."""
    SUBSONIC = "Subsonic"           # M < 0.8
    TRANSONIC = "Transonic"         # 0.8 <= M <= 1.2
    SUPERSONIC = "Supersonic"       # 1.2 < M < 5.0
    HYPERSONIC = "Hypersonic"       # M >= 5.0


def get_flow_regime(mach: float) -> FlowRegime:
    """Classify flow regime based on Mach number."""
    if mach < 0.8:
        return FlowRegime.SUBSONIC
    elif mach <= 1.2:
        return FlowRegime.TRANSONIC
    elif mach < 5.0:
        return FlowRegime.SUPERSONIC
    else:
        return FlowRegime.HYPERSONIC


@dataclass
class WingGeometry:
    """Aircraft wing geometry parameters."""
    span: float                      # Wingspan b (m)
    root_chord: float                # Root chord c_r (m)
    tip_chord: float                 # Tip chord c_t (m)
    sweep_leading_edge: float = 0.0  # Leading edge sweep angle (rad)

    @property
    def taper_ratio(self) -> float:
        """Taper ratio λ = c_tip / c_root"""
        return self.tip_chord / self.root_chord

    @property
    def mean_aerodynamic_chord(self) -> float:
        """Mean aerodynamic chord (MAC)"""
        lam = self.taper_ratio
        return (2/3) * self.root_chord * (1 + lam + lam**2) / (1 + lam)

    @property
    def wing_area(self) -> float:
        """Wing planform area S"""
        return 0.5 * self.span * (self.root_chord + self.tip_chord)

    @property
    def aspect_ratio(self) -> float:
        """Aspect ratio AR = b²/S"""
        return self.span**2 / self.wing_area

    @property
    def sweep_quarter_chord(self) -> float:
        """Quarter-chord sweep angle (approximate conversion from LE sweep)"""
        lam = self.taper_ratio
        ar = self.aspect_ratio
        # Conversion formula from LE sweep to c/4 sweep
        tan_le = np.tan(self.sweep_leading_edge)
        tan_c4 = tan_le - (1 - lam) / (ar * (1 + lam))
        return np.arctan(tan_c4)


def mach_cone_angle(mach: float) -> float:
    """
    Calculate Mach cone half-angle (μ) for supersonic flow.

    μ = arcsin(1/M)

    Parameters:
        mach: Freestream Mach number (must be > 1)

    Returns:
        mu: Mach cone half-angle (rad)
    """
    if mach <= 1.0:
        return np.pi / 2  # 90° for subsonic (no cone)
    return np.arcsin(1.0 / mach)


def is_wing_subsonic_le(mach: float, sweep_le: float) -> bool:
    """
    Check if wing has subsonic or supersonic leading edge.

    Subsonic LE: Λ_LE > μ (sweep greater than Mach angle)
    Supersonic LE: Λ_LE < μ (sweep less than Mach angle)

    Parameters:
        mach: Freestream Mach number
        sweep_le: Leading edge sweep angle (rad)

    Returns:
        True if subsonic leading edge, False if supersonic
    """
    if mach <= 1.0:
        return True
    mu = mach_cone_angle(mach)
    return sweep_le > mu


def compressibility_factor(mach: float) -> float:
    """
    Calculate compressibility factor β for different flow regimes.

    Subsonic:   β = √(1 - M²)     (Prandtl-Glauert)
    Supersonic: β = √(M² - 1)     (Ackeret)

    Parameters:
        mach: Freestream Mach number

    Returns:
        beta: Compressibility factor
    """
    if mach < 1.0:
        return np.sqrt(1.0 - mach**2)
    else:
        return np.sqrt(mach**2 - 1.0)


# =============================================================================
# LIFT CURVE SLOPE CALCULATIONS
# =============================================================================

def lift_curve_slope_2d_subsonic(mach: float = 0.0) -> float:
    """
    2D airfoil lift curve slope for subsonic flow.

    C_lα = 2π / β  (Prandtl-Glauert corrected thin airfoil theory)

    Parameters:
        mach: Freestream Mach number

    Returns:
        cl_alpha_2d: 2D lift curve slope (1/rad)
    """
    beta = compressibility_factor(mach) if mach < 1.0 else 1.0
    return 2 * np.pi / beta


def lift_curve_slope_2d_supersonic(mach: float) -> float:
    """
    2D airfoil lift curve slope for supersonic flow (Ackeret theory).

    C_lα = 4 / √(M² - 1) = 4 / β

    Parameters:
        mach: Freestream Mach number (must be > 1)

    Returns:
        cl_alpha_2d: 2D lift curve slope (1/rad)
    """
    if mach <= 1.0:
        raise ValueError("Mach number must be > 1 for supersonic theory")
    beta = compressibility_factor(mach)
    return 4.0 / beta


def lift_curve_slope_3d_subsonic(aspect_ratio: float, sweep: float = 0.0,
                                  mach: float = 0.0) -> float:
    """
    3D wing lift curve slope for subsonic flow (Helmbold/Polhamus).

    Parameters:
        aspect_ratio: Wing aspect ratio AR
        sweep: Quarter-chord sweep angle (rad)
        mach: Freestream Mach number

    Returns:
        cl_alpha_3d: 3D wing lift curve slope (1/rad)
    """
    cl_alpha_2d = lift_curve_slope_2d_subsonic(mach)
    beta = compressibility_factor(mach) if mach < 1.0 else 1.0
    cos_sweep = np.cos(sweep)

    # Effective aspect ratio with sweep
    ar_eff = aspect_ratio * cos_sweep

    # Helmbold equation
    cl_alpha_3d = (cl_alpha_2d * ar_eff) / (
        2 + np.sqrt(4 + (ar_eff * beta / cl_alpha_2d * cos_sweep)**2)
    )

    return cl_alpha_3d


def lift_curve_slope_3d_supersonic(wing: WingGeometry, mach: float) -> float:
    """
    3D wing lift curve slope for supersonic flow.

    For supersonic leading edge (unswept or low sweep):
        C_Lα = 4 / √(M² - 1)  (2D Ackeret, no 3D correction needed for high AR)

    For subsonic leading edge (highly swept):
        Uses swept wing supersonic theory with reduced effectiveness

    Parameters:
        wing: WingGeometry object
        mach: Freestream Mach number (must be > 1)

    Returns:
        cl_alpha_3d: 3D wing lift curve slope (1/rad)
    """
    if mach <= 1.0:
        raise ValueError("Mach number must be > 1 for supersonic calculation")

    beta = compressibility_factor(mach)
    ar = wing.aspect_ratio
    sweep_le = wing.sweep_leading_edge

    # Check leading edge type
    if is_wing_subsonic_le(mach, sweep_le):
        # Subsonic leading edge - use modified theory
        # Jones' slender wing theory / subsonic LE correction
        m = beta * np.tan(sweep_le)  # sweep parameter

        if m > 1:
            # Highly swept - subsonic LE dominates
            cos_sweep = np.cos(sweep_le)
            # Effective Mach number normal to LE
            mach_n = mach * cos_sweep
            if mach_n < 1.0:
                # Treat as subsonic in normal direction
                beta_n = np.sqrt(1 - mach_n**2)
                cl_alpha_3d = (2 * np.pi * ar) / (
                    2 + np.sqrt(ar**2 * beta_n**2 + 4)
                )
            else:
                beta_n = np.sqrt(mach_n**2 - 1)
                cl_alpha_3d = 4 / beta_n * (ar / (ar + 2))
        else:
            cl_alpha_3d = 4 / beta
    else:
        # Supersonic leading edge
        # For high aspect ratio: approaches 2D value
        # For finite AR: apply correction
        cl_alpha_2d = 4.0 / beta

        # Finite wing correction for supersonic flow
        # Based on linear supersonic wing theory
        if ar * beta > 4:
            # High AR approximation
            cl_alpha_3d = cl_alpha_2d * (1 - 1/(ar * beta))
        else:
            # Low AR correction
            cl_alpha_3d = (np.pi * ar / 2) * (1 / (1 + np.sqrt(1 + (beta * ar / 4)**2)))

    return cl_alpha_3d


def lift_curve_slope_3d(wing: WingGeometry, mach: float) -> float:
    """
    Unified 3D lift curve slope calculation for all flow regimes.

    Parameters:
        wing: WingGeometry object
        mach: Freestream Mach number

    Returns:
        cl_alpha_3d: 3D wing lift curve slope (1/rad)
    """
    regime = get_flow_regime(mach)

    if regime == FlowRegime.SUBSONIC:
        return lift_curve_slope_3d_subsonic(
            wing.aspect_ratio, wing.sweep_quarter_chord, mach
        )

    elif regime == FlowRegime.TRANSONIC:
        # Linear interpolation between subsonic and supersonic
        # (simplified - real transonic is highly nonlinear)
        cl_sub = lift_curve_slope_3d_subsonic(
            wing.aspect_ratio, wing.sweep_quarter_chord, 0.8
        )
        cl_sup = lift_curve_slope_3d_supersonic(wing, 1.2)

        # Interpolation factor
        t = (mach - 0.8) / 0.4
        return cl_sub * (1 - t) + cl_sup * t

    else:  # SUPERSONIC or HYPERSONIC
        return lift_curve_slope_3d_supersonic(wing, mach)


# =============================================================================
# ROLL DAMPING DERIVATIVE CALCULATIONS
# =============================================================================

def calculate_clp_subsonic(wing: WingGeometry, mach: float = 0.0) -> float:
    """
    Roll damping derivative for subsonic flow.

    C_l_p = -(C_Lα/12) * (1 + 3λ) / (1 + λ) * K_sweep

    Parameters:
        wing: WingGeometry object
        mach: Freestream Mach number (< 1)

    Returns:
        clp: Roll damping derivative (1/rad)
    """
    cl_alpha = lift_curve_slope_3d_subsonic(
        wing.aspect_ratio, wing.sweep_quarter_chord, mach
    )
    lam = wing.taper_ratio
    sweep = wing.sweep_quarter_chord

    # Taper ratio factor
    taper_factor = (1 + 3*lam) / (1 + lam)

    # Sweep correction factor
    sweep_factor = np.cos(sweep)**2

    return -(cl_alpha / 12) * taper_factor * sweep_factor


def calculate_clp_supersonic_analytical(wing: WingGeometry, mach: float) -> float:
    """
    Roll damping derivative for supersonic flow using analytical method.

    For supersonic leading edge (Ackeret-based):
        C_l_p = -4/(3β) * (AR/(AR + 4/β))

    For subsonic leading edge (swept wing):
        Modified formula accounting for sweep effects

    Parameters:
        wing: WingGeometry object
        mach: Freestream Mach number (> 1)

    Returns:
        clp: Roll damping derivative (1/rad)
    """
    if mach <= 1.0:
        raise ValueError("Mach number must be > 1 for supersonic calculation")

    beta = compressibility_factor(mach)
    ar = wing.aspect_ratio
    lam = wing.taper_ratio
    sweep_le = wing.sweep_leading_edge

    # Taper ratio factor (modified for supersonic)
    taper_factor = (1 + 2*lam) / (1 + lam)

    if is_wing_subsonic_le(mach, sweep_le):
        # Subsonic leading edge - reduced effectiveness
        cos_sweep = np.cos(sweep_le)
        mach_n = mach * cos_sweep

        if mach_n < 1.0:
            # Subsonic normal Mach - use modified subsonic formula
            beta_n = np.sqrt(1 - mach_n**2)
            clp = -(2 * np.pi / (3 * (ar + 2/beta_n))) * taper_factor
        else:
            # Supersonic normal Mach
            beta_n = np.sqrt(mach_n**2 - 1)
            clp = -(4 / (3 * beta_n)) * (ar / (ar + 4/beta_n)) * taper_factor

        # Apply sweep reduction
        clp *= cos_sweep**2

    else:
        # Supersonic leading edge - full supersonic theory
        clp = -(4 / (3 * beta)) * (ar / (ar + 4/beta)) * taper_factor

    return clp


def calculate_clp_supersonic_strip(wing: WingGeometry, mach: float,
                                    n_strips: int = 100) -> float:
    """
    Roll damping derivative for supersonic flow using strip theory.

    Integrates using local supersonic airfoil theory along span.

    Parameters:
        wing: WingGeometry object
        mach: Freestream Mach number (> 1)
        n_strips: Number of integration strips

    Returns:
        clp: Roll damping derivative (1/rad)
    """
    if mach <= 1.0:
        raise ValueError("Mach number must be > 1 for supersonic calculation")

    b = wing.span
    S = wing.wing_area
    c_r = wing.root_chord
    lam = wing.taper_ratio
    beta = compressibility_factor(mach)
    sweep_le = wing.sweep_leading_edge

    # Spanwise stations
    y = np.linspace(0.01, b/2, n_strips)  # Avoid y=0 singularity

    # Local chord distribution
    c_y = c_r * (1 - (1 - lam) * 2 * y / b)

    # Local lift curve slope (supersonic 2D theory)
    # Modified for 3D tip effects
    tip_loss = 1 - (2 * y / b)**2  # Approximate tip loss factor
    cl_alpha_local = (4.0 / beta) * tip_loss

    # Sweep effect on local flow
    if sweep_le > 0:
        cos_sweep = np.cos(sweep_le)
        mach_n = mach * cos_sweep
        if mach_n > 1:
            beta_n = np.sqrt(mach_n**2 - 1)
            cl_alpha_local = (4.0 / beta_n) * tip_loss * cos_sweep

    # Integration
    integrand = c_y * cl_alpha_local * y**2
    clp = -(4 / (S * b**2)) * np.trapezoid(integrand, y)

    return clp


def calculate_clp(wing: WingGeometry, mach: float) -> Tuple[float, FlowRegime]:
    """
    Unified roll damping derivative calculation for all flow regimes.

    Parameters:
        wing: WingGeometry object
        mach: Freestream Mach number

    Returns:
        clp: Roll damping derivative (1/rad)
        regime: Flow regime classification
    """
    regime = get_flow_regime(mach)

    if regime == FlowRegime.SUBSONIC:
        clp = calculate_clp_subsonic(wing, mach)

    elif regime == FlowRegime.TRANSONIC:
        # Interpolation through transonic region
        clp_sub = calculate_clp_subsonic(wing, 0.8)
        clp_sup = calculate_clp_supersonic_analytical(wing, 1.2)
        t = (mach - 0.8) / 0.4
        # Nonlinear interpolation (transonic dip)
        t_smooth = 3*t**2 - 2*t**3  # Smoothstep
        clp = clp_sub * (1 - t_smooth) + clp_sup * t_smooth

    elif regime == FlowRegime.SUPERSONIC:
        clp = calculate_clp_supersonic_analytical(wing, mach)

    else:  # HYPERSONIC
        # Use supersonic with Mach independence approximation
        clp = calculate_clp_supersonic_analytical(wing, 5.0)
        # Newtonian correction for very high Mach
        clp *= (5.0 / mach)**0.5

    return clp, regime


def dimensional_roll_damping(clp: float, rho: float, V: float,
                              wing: WingGeometry) -> float:
    """
    Convert non-dimensional C_l_p to dimensional roll damping derivative L_p.

    L_p = ∂L/∂p = C_l_p * (ρVSb²)/(4)

    Parameters:
        clp: Non-dimensional roll damping derivative (1/rad)
        rho: Air density (kg/m³)
        V: True airspeed (m/s)
        wing: WingGeometry object

    Returns:
        Lp: Dimensional roll damping (N·m·s/rad)
    """
    S = wing.wing_area
    b = wing.span
    return clp * (rho * V * S * b**2) / 4


def standard_atmosphere(altitude: float) -> Tuple[float, float, float]:
    """
    ISA standard atmosphere properties.

    Parameters:
        altitude: Geometric altitude (m)

    Returns:
        temperature: Temperature (K)
        pressure: Pressure (Pa)
        density: Density (kg/m³)
    """
    if altitude < 11000:
        # Troposphere
        T = 288.15 - 0.0065 * altitude
        p = 101325 * (T / 288.15)**5.2561
    else:
        # Stratosphere (up to 20km)
        T = 216.65
        p = 22632 * np.exp(-9.81 * (altitude - 11000) / (287 * T))

    rho = p / (287 * T)
    return T, p, rho


def main():
    """Example calculations for subsonic and supersonic flight."""

    print("=" * 70)
    print("Roll Damping Derivative (C_l_p) Calculator")
    print("SUBSONIC and SUPERSONIC Regimes")
    print("=" * 70)

    # Example 1: Subsonic - General Aviation (Cessna-like)
    print("\n" + "=" * 70)
    print("EXAMPLE 1: SUBSONIC - General Aviation Aircraft")
    print("=" * 70)

    wing_subsonic = WingGeometry(
        span=11.0,
        root_chord=1.63,
        tip_chord=1.12,
        sweep_leading_edge=np.radians(2)
    )

    mach_sub = 0.2
    altitude_sub = 3000
    _, _, rho_sub = standard_atmosphere(altitude_sub)

    print(f"\nWing Geometry:")
    print(f"  Wingspan:         {wing_subsonic.span:.2f} m")
    print(f"  Taper ratio:      {wing_subsonic.taper_ratio:.3f}")
    print(f"  Aspect ratio:     {wing_subsonic.aspect_ratio:.2f}")
    print(f"  LE Sweep:         {np.degrees(wing_subsonic.sweep_leading_edge):.1f}°")

    print(f"\nFlight Conditions:")
    print(f"  Mach number:      {mach_sub}")
    print(f"  Altitude:         {altitude_sub} m")

    clp_sub, regime_sub = calculate_clp(wing_subsonic, mach_sub)
    cl_alpha_sub = lift_curve_slope_3d(wing_subsonic, mach_sub)

    print(f"\nResults ({regime_sub.value}):")
    print(f"  C_Lα (3D):        {cl_alpha_sub:.4f} 1/rad")
    print(f"  C_l_p:            {clp_sub:.4f} 1/rad")

    # Example 2: Supersonic - Fighter Aircraft (F-16-like)
    print("\n" + "=" * 70)
    print("EXAMPLE 2: SUPERSONIC - Fighter Aircraft")
    print("=" * 70)

    wing_supersonic = WingGeometry(
        span=9.45,
        root_chord=5.0,
        tip_chord=1.5,
        sweep_leading_edge=np.radians(40)
    )

    mach_sup = 1.6
    altitude_sup = 12000
    _, _, rho_sup = standard_atmosphere(altitude_sup)

    print(f"\nWing Geometry:")
    print(f"  Wingspan:         {wing_supersonic.span:.2f} m")
    print(f"  Taper ratio:      {wing_supersonic.taper_ratio:.3f}")
    print(f"  Aspect ratio:     {wing_supersonic.aspect_ratio:.2f}")
    print(f"  LE Sweep:         {np.degrees(wing_supersonic.sweep_leading_edge):.1f}°")

    mu = mach_cone_angle(mach_sup)
    le_type = "Subsonic LE" if is_wing_subsonic_le(mach_sup, wing_supersonic.sweep_leading_edge) else "Supersonic LE"

    print(f"\nFlight Conditions:")
    print(f"  Mach number:      {mach_sup}")
    print(f"  Altitude:         {altitude_sup} m")
    print(f"  Mach angle (μ):   {np.degrees(mu):.1f}°")
    print(f"  Leading edge:     {le_type}")

    clp_sup, regime_sup = calculate_clp(wing_supersonic, mach_sup)
    cl_alpha_sup = lift_curve_slope_3d(wing_supersonic, mach_sup)
    beta_sup = compressibility_factor(mach_sup)

    print(f"\nResults ({regime_sup.value}):")
    print(f"  β = √(M²-1):      {beta_sup:.4f}")
    print(f"  C_Lα (3D):        {cl_alpha_sup:.4f} 1/rad")
    print(f"  C_l_p:            {clp_sup:.4f} 1/rad")

    # Strip theory comparison
    clp_strip = calculate_clp_supersonic_strip(wing_supersonic, mach_sup)
    print(f"  C_l_p (strip):    {clp_strip:.4f} 1/rad")

    # Example 3: Mach sweep
    print("\n" + "=" * 70)
    print("EXAMPLE 3: C_l_p vs MACH NUMBER")
    print("=" * 70)

    mach_range = [0.3, 0.6, 0.8, 0.9, 1.0, 1.1, 1.2, 1.5, 2.0, 2.5, 3.0]

    print(f"\nWing: AR={wing_supersonic.aspect_ratio:.2f}, "
          f"λ={wing_supersonic.taper_ratio:.2f}, "
          f"Λ_LE={np.degrees(wing_supersonic.sweep_leading_edge):.0f}°")
    print("-" * 50)
    print(f"{'Mach':>6} | {'Regime':^12} | {'C_Lα':>8} | {'C_l_p':>8}")
    print("-" * 50)

    for M in mach_range:
        clp, regime = calculate_clp(wing_supersonic, M)
        cl_a = lift_curve_slope_3d(wing_supersonic, M)
        print(f"{M:>6.2f} | {regime.value:^12} | {cl_a:>8.4f} | {clp:>8.4f}")

    print("-" * 50)
    print("\nNotes:")
    print("  - Negative C_l_p indicates roll damping (stable)")
    print("  - |C_l_p| generally decreases with increasing Mach")
    print("  - Transonic region uses interpolation (simplified)")
    print("  - Subsonic LE provides better supersonic damping")


if __name__ == "__main__":
    main()
