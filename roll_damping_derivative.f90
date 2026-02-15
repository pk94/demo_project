!==============================================================================
! Roll Damping Derivative Calculator (C_l_p)
!==============================================================================
! Calculates the derivative of rolling moment coefficient with respect to
! non-dimensional roll angular velocity for aircraft wings.
!
! Supports: SUBSONIC, TRANSONIC, and SUPERSONIC flight regimes.
!
! In flight dynamics notation:
!     C_l_p = dC_l / dp_hat
!
! where p_hat = pb/(2V) is the non-dimensional roll rate
!     p = roll angular velocity (rad/s)
!     b = wingspan (m)
!     V = true airspeed (m/s)
!
! Compile: gfortran -o roll_damping roll_damping_derivative.f90
! Run:     ./roll_damping
!==============================================================================

module constants
    implicit none

    ! Double precision kind
    integer, parameter :: dp = selected_real_kind(15, 307)

    ! Mathematical constants
    real(dp), parameter :: PI = 3.141592653589793_dp
    real(dp), parameter :: DEG2RAD = PI / 180.0_dp
    real(dp), parameter :: RAD2DEG = 180.0_dp / PI

    ! Flow regime identifiers
    integer, parameter :: REGIME_SUBSONIC   = 1
    integer, parameter :: REGIME_TRANSONIC  = 2
    integer, parameter :: REGIME_SUPERSONIC = 3
    integer, parameter :: REGIME_HYPERSONIC = 4

end module constants


module wing_geometry_mod
    use constants
    implicit none

    ! Wing geometry derived type
    type :: WingGeometry
        real(dp) :: span              ! Wingspan b (m)
        real(dp) :: root_chord        ! Root chord c_r (m)
        real(dp) :: tip_chord         ! Tip chord c_t (m)
        real(dp) :: sweep_le          ! Leading edge sweep angle (rad)
    end type WingGeometry

contains

    !--------------------------------------------------------------------------
    ! Calculate taper ratio: lambda = c_tip / c_root
    !--------------------------------------------------------------------------
    pure function taper_ratio(wing) result(lambda)
        type(WingGeometry), intent(in) :: wing
        real(dp) :: lambda

        lambda = wing%tip_chord / wing%root_chord
    end function taper_ratio

    !--------------------------------------------------------------------------
    ! Calculate wing planform area: S = (b/2) * (c_root + c_tip)
    !--------------------------------------------------------------------------
    pure function wing_area(wing) result(S)
        type(WingGeometry), intent(in) :: wing
        real(dp) :: S

        S = 0.5_dp * wing%span * (wing%root_chord + wing%tip_chord)
    end function wing_area

    !--------------------------------------------------------------------------
    ! Calculate aspect ratio: AR = b^2 / S
    !--------------------------------------------------------------------------
    pure function aspect_ratio(wing) result(AR)
        type(WingGeometry), intent(in) :: wing
        real(dp) :: AR, S

        S = wing_area(wing)
        AR = wing%span**2 / S
    end function aspect_ratio

    !--------------------------------------------------------------------------
    ! Calculate mean aerodynamic chord (MAC)
    !--------------------------------------------------------------------------
    pure function mean_aero_chord(wing) result(MAC)
        type(WingGeometry), intent(in) :: wing
        real(dp) :: MAC, lam

        lam = taper_ratio(wing)
        MAC = (2.0_dp/3.0_dp) * wing%root_chord * &
              (1.0_dp + lam + lam**2) / (1.0_dp + lam)
    end function mean_aero_chord

    !--------------------------------------------------------------------------
    ! Calculate quarter-chord sweep from leading edge sweep
    !--------------------------------------------------------------------------
    pure function sweep_quarter_chord(wing) result(sweep_c4)
        type(WingGeometry), intent(in) :: wing
        real(dp) :: sweep_c4, lam, AR, tan_le, tan_c4

        lam = taper_ratio(wing)
        AR = aspect_ratio(wing)
        tan_le = tan(wing%sweep_le)
        tan_c4 = tan_le - (1.0_dp - lam) / (AR * (1.0_dp + lam))
        sweep_c4 = atan(tan_c4)
    end function sweep_quarter_chord

end module wing_geometry_mod


module flow_regime_mod
    use constants
    implicit none

contains

    !--------------------------------------------------------------------------
    ! Classify flow regime based on Mach number
    !--------------------------------------------------------------------------
    pure function get_flow_regime(mach) result(regime)
        real(dp), intent(in) :: mach
        integer :: regime

        if (mach < 0.8_dp) then
            regime = REGIME_SUBSONIC
        else if (mach <= 1.2_dp) then
            regime = REGIME_TRANSONIC
        else if (mach < 5.0_dp) then
            regime = REGIME_SUPERSONIC
        else
            regime = REGIME_HYPERSONIC
        end if
    end function get_flow_regime

    !--------------------------------------------------------------------------
    ! Get flow regime name as string
    !--------------------------------------------------------------------------
    function regime_name(regime) result(name)
        integer, intent(in) :: regime
        character(len=12) :: name

        select case (regime)
            case (REGIME_SUBSONIC)
                name = "Subsonic"
            case (REGIME_TRANSONIC)
                name = "Transonic"
            case (REGIME_SUPERSONIC)
                name = "Supersonic"
            case (REGIME_HYPERSONIC)
                name = "Hypersonic"
            case default
                name = "Unknown"
        end select
    end function regime_name

    !--------------------------------------------------------------------------
    ! Calculate Mach cone half-angle: mu = arcsin(1/M)
    !--------------------------------------------------------------------------
    pure function mach_cone_angle(mach) result(mu)
        real(dp), intent(in) :: mach
        real(dp) :: mu

        if (mach <= 1.0_dp) then
            mu = PI / 2.0_dp  ! 90 degrees for subsonic
        else
            mu = asin(1.0_dp / mach)
        end if
    end function mach_cone_angle

    !--------------------------------------------------------------------------
    ! Check if wing has subsonic leading edge
    ! Subsonic LE: sweep_LE > mu (Mach angle)
    !--------------------------------------------------------------------------
    pure function is_subsonic_le(mach, sweep_le) result(is_sub)
        real(dp), intent(in) :: mach, sweep_le
        logical :: is_sub
        real(dp) :: mu

        if (mach <= 1.0_dp) then
            is_sub = .true.
        else
            mu = mach_cone_angle(mach)
            is_sub = (sweep_le > mu)
        end if
    end function is_subsonic_le

    !--------------------------------------------------------------------------
    ! Calculate compressibility factor beta
    ! Subsonic:   beta = sqrt(1 - M^2)
    ! Supersonic: beta = sqrt(M^2 - 1)
    !--------------------------------------------------------------------------
    pure function compressibility_factor(mach) result(beta)
        real(dp), intent(in) :: mach
        real(dp) :: beta

        if (mach < 1.0_dp) then
            beta = sqrt(1.0_dp - mach**2)
        else
            beta = sqrt(mach**2 - 1.0_dp)
        end if
    end function compressibility_factor

end module flow_regime_mod


module lift_curve_slope_mod
    use constants
    use wing_geometry_mod
    use flow_regime_mod
    implicit none

contains

    !--------------------------------------------------------------------------
    ! 2D airfoil lift curve slope - Subsonic
    ! C_l_alpha = 2*pi / beta (Prandtl-Glauert corrected)
    !--------------------------------------------------------------------------
    pure function cl_alpha_2d_subsonic(mach) result(cl_alpha)
        real(dp), intent(in) :: mach
        real(dp) :: cl_alpha, beta

        if (mach < 1.0_dp) then
            beta = compressibility_factor(mach)
        else
            beta = 1.0_dp
        end if

        cl_alpha = 2.0_dp * PI / beta
    end function cl_alpha_2d_subsonic

    !--------------------------------------------------------------------------
    ! 2D airfoil lift curve slope - Supersonic (Ackeret theory)
    ! C_l_alpha = 4 / beta = 4 / sqrt(M^2 - 1)
    !--------------------------------------------------------------------------
    pure function cl_alpha_2d_supersonic(mach) result(cl_alpha)
        real(dp), intent(in) :: mach
        real(dp) :: cl_alpha, beta

        beta = compressibility_factor(mach)
        cl_alpha = 4.0_dp / beta
    end function cl_alpha_2d_supersonic

    !--------------------------------------------------------------------------
    ! 3D wing lift curve slope - Subsonic (Helmbold equation)
    !--------------------------------------------------------------------------
    pure function cl_alpha_3d_subsonic(AR, sweep, mach) result(cl_alpha)
        real(dp), intent(in) :: AR, sweep, mach
        real(dp) :: cl_alpha, cl_alpha_2d, beta, cos_sweep, AR_eff

        cl_alpha_2d = cl_alpha_2d_subsonic(mach)

        if (mach < 1.0_dp) then
            beta = compressibility_factor(mach)
        else
            beta = 1.0_dp
        end if

        cos_sweep = cos(sweep)
        AR_eff = AR * cos_sweep

        ! Helmbold equation
        cl_alpha = (cl_alpha_2d * AR_eff) / &
                   (2.0_dp + sqrt(4.0_dp + (AR_eff * beta / cl_alpha_2d * cos_sweep)**2))
    end function cl_alpha_3d_subsonic

    !--------------------------------------------------------------------------
    ! 3D wing lift curve slope - Supersonic
    !--------------------------------------------------------------------------
    pure function cl_alpha_3d_supersonic(wing, mach) result(cl_alpha)
        type(WingGeometry), intent(in) :: wing
        real(dp), intent(in) :: mach
        real(dp) :: cl_alpha
        real(dp) :: beta, AR, sweep_le, m_param
        real(dp) :: cos_sweep, mach_n, beta_n, cl_alpha_2d

        beta = compressibility_factor(mach)
        AR = aspect_ratio(wing)
        sweep_le = wing%sweep_le

        if (is_subsonic_le(mach, sweep_le)) then
            ! Subsonic leading edge
            m_param = beta * tan(sweep_le)

            if (m_param > 1.0_dp) then
                ! Highly swept - subsonic LE dominates
                cos_sweep = cos(sweep_le)
                mach_n = mach * cos_sweep

                if (mach_n < 1.0_dp) then
                    ! Subsonic in normal direction
                    beta_n = sqrt(1.0_dp - mach_n**2)
                    cl_alpha = (2.0_dp * PI * AR) / &
                               (2.0_dp + sqrt(AR**2 * beta_n**2 + 4.0_dp))
                else
                    ! Supersonic normal Mach
                    beta_n = sqrt(mach_n**2 - 1.0_dp)
                    cl_alpha = 4.0_dp / beta_n * (AR / (AR + 2.0_dp))
                end if
            else
                cl_alpha = 4.0_dp / beta
            end if
        else
            ! Supersonic leading edge
            cl_alpha_2d = 4.0_dp / beta

            if (AR * beta > 4.0_dp) then
                ! High AR approximation
                cl_alpha = cl_alpha_2d * (1.0_dp - 1.0_dp/(AR * beta))
            else
                ! Low AR correction
                cl_alpha = (PI * AR / 2.0_dp) * &
                           (1.0_dp / (1.0_dp + sqrt(1.0_dp + (beta * AR / 4.0_dp)**2)))
            end if
        end if
    end function cl_alpha_3d_supersonic

    !--------------------------------------------------------------------------
    ! Unified 3D lift curve slope for all regimes
    !--------------------------------------------------------------------------
    pure function cl_alpha_3d(wing, mach) result(cl_alpha)
        type(WingGeometry), intent(in) :: wing
        real(dp), intent(in) :: mach
        real(dp) :: cl_alpha
        integer :: regime
        real(dp) :: cl_sub, cl_sup, t, AR, sweep

        regime = get_flow_regime(mach)
        AR = aspect_ratio(wing)
        sweep = sweep_quarter_chord(wing)

        select case (regime)
            case (REGIME_SUBSONIC)
                cl_alpha = cl_alpha_3d_subsonic(AR, sweep, mach)

            case (REGIME_TRANSONIC)
                ! Linear interpolation
                cl_sub = cl_alpha_3d_subsonic(AR, sweep, 0.8_dp)
                cl_sup = cl_alpha_3d_supersonic(wing, 1.2_dp)
                t = (mach - 0.8_dp) / 0.4_dp
                cl_alpha = cl_sub * (1.0_dp - t) + cl_sup * t

            case default  ! SUPERSONIC or HYPERSONIC
                cl_alpha = cl_alpha_3d_supersonic(wing, mach)
        end select
    end function cl_alpha_3d

end module lift_curve_slope_mod


module roll_damping_mod
    use constants
    use wing_geometry_mod
    use flow_regime_mod
    use lift_curve_slope_mod
    implicit none

contains

    !--------------------------------------------------------------------------
    ! Roll damping derivative - Subsonic
    ! C_l_p = -(C_L_alpha/12) * (1 + 3*lambda)/(1 + lambda) * cos^2(sweep)
    !--------------------------------------------------------------------------
    pure function clp_subsonic(wing, mach) result(clp)
        type(WingGeometry), intent(in) :: wing
        real(dp), intent(in) :: mach
        real(dp) :: clp
        real(dp) :: cl_alpha, lam, sweep, taper_factor, sweep_factor

        cl_alpha = cl_alpha_3d_subsonic(aspect_ratio(wing), &
                                        sweep_quarter_chord(wing), mach)
        lam = taper_ratio(wing)
        sweep = sweep_quarter_chord(wing)

        taper_factor = (1.0_dp + 3.0_dp*lam) / (1.0_dp + lam)
        sweep_factor = cos(sweep)**2

        clp = -(cl_alpha / 12.0_dp) * taper_factor * sweep_factor
    end function clp_subsonic

    !--------------------------------------------------------------------------
    ! Roll damping derivative - Supersonic (Analytical)
    !--------------------------------------------------------------------------
    pure function clp_supersonic_analytical(wing, mach) result(clp)
        type(WingGeometry), intent(in) :: wing
        real(dp), intent(in) :: mach
        real(dp) :: clp
        real(dp) :: beta, AR, lam, sweep_le, taper_factor
        real(dp) :: cos_sweep, mach_n, beta_n

        beta = compressibility_factor(mach)
        AR = aspect_ratio(wing)
        lam = taper_ratio(wing)
        sweep_le = wing%sweep_le

        ! Modified taper factor for supersonic
        taper_factor = (1.0_dp + 2.0_dp*lam) / (1.0_dp + lam)

        if (is_subsonic_le(mach, sweep_le)) then
            ! Subsonic leading edge
            cos_sweep = cos(sweep_le)
            mach_n = mach * cos_sweep

            if (mach_n < 1.0_dp) then
                ! Subsonic normal Mach
                beta_n = sqrt(1.0_dp - mach_n**2)
                clp = -(2.0_dp * PI / (3.0_dp * (AR + 2.0_dp/beta_n))) * taper_factor
            else
                ! Supersonic normal Mach
                beta_n = sqrt(mach_n**2 - 1.0_dp)
                clp = -(4.0_dp / (3.0_dp * beta_n)) * &
                      (AR / (AR + 4.0_dp/beta_n)) * taper_factor
            end if

            ! Apply sweep reduction
            clp = clp * cos_sweep**2
        else
            ! Supersonic leading edge
            clp = -(4.0_dp / (3.0_dp * beta)) * &
                  (AR / (AR + 4.0_dp/beta)) * taper_factor
        end if
    end function clp_supersonic_analytical

    !--------------------------------------------------------------------------
    ! Roll damping derivative - Supersonic (Strip Theory)
    !--------------------------------------------------------------------------
    pure function clp_supersonic_strip(wing, mach, n_strips) result(clp)
        type(WingGeometry), intent(in) :: wing
        real(dp), intent(in) :: mach
        integer, intent(in) :: n_strips
        real(dp) :: clp

        real(dp) :: b, S, c_r, lam, beta, sweep_le
        real(dp) :: y, dy, c_y, tip_loss, cl_alpha_local
        real(dp) :: cos_sweep, mach_n, beta_n
        real(dp) :: integral
        integer :: i

        b = wing%span
        S = wing_area(wing)
        c_r = wing%root_chord
        lam = taper_ratio(wing)
        beta = compressibility_factor(mach)
        sweep_le = wing%sweep_le

        dy = (b / 2.0_dp - 0.01_dp) / real(n_strips - 1, dp)
        integral = 0.0_dp

        do i = 1, n_strips
            y = 0.01_dp + real(i - 1, dp) * dy

            ! Local chord
            c_y = c_r * (1.0_dp - (1.0_dp - lam) * 2.0_dp * y / b)

            ! Tip loss factor
            tip_loss = 1.0_dp - (2.0_dp * y / b)**2

            ! Local lift curve slope
            cl_alpha_local = (4.0_dp / beta) * tip_loss

            ! Sweep effect
            if (sweep_le > 0.0_dp) then
                cos_sweep = cos(sweep_le)
                mach_n = mach * cos_sweep
                if (mach_n > 1.0_dp) then
                    beta_n = sqrt(mach_n**2 - 1.0_dp)
                    cl_alpha_local = (4.0_dp / beta_n) * tip_loss * cos_sweep
                end if
            end if

            ! Trapezoidal integration
            if (i == 1 .or. i == n_strips) then
                integral = integral + 0.5_dp * c_y * cl_alpha_local * y**2 * dy
            else
                integral = integral + c_y * cl_alpha_local * y**2 * dy
            end if
        end do

        clp = -(4.0_dp / (S * b**2)) * integral
    end function clp_supersonic_strip

    !--------------------------------------------------------------------------
    ! Unified roll damping derivative for all regimes
    !--------------------------------------------------------------------------
    pure subroutine calculate_clp(wing, mach, clp, regime)
        type(WingGeometry), intent(in) :: wing
        real(dp), intent(in) :: mach
        real(dp), intent(out) :: clp
        integer, intent(out) :: regime

        real(dp) :: clp_sub, clp_sup, t, t_smooth

        regime = get_flow_regime(mach)

        select case (regime)
            case (REGIME_SUBSONIC)
                clp = clp_subsonic(wing, mach)

            case (REGIME_TRANSONIC)
                ! Interpolation with smoothstep
                clp_sub = clp_subsonic(wing, 0.8_dp)
                clp_sup = clp_supersonic_analytical(wing, 1.2_dp)
                t = (mach - 0.8_dp) / 0.4_dp
                t_smooth = 3.0_dp * t**2 - 2.0_dp * t**3
                clp = clp_sub * (1.0_dp - t_smooth) + clp_sup * t_smooth

            case (REGIME_SUPERSONIC)
                clp = clp_supersonic_analytical(wing, mach)

            case (REGIME_HYPERSONIC)
                ! Mach independence with Newtonian correction
                clp = clp_supersonic_analytical(wing, 5.0_dp)
                clp = clp * sqrt(5.0_dp / mach)
        end select
    end subroutine calculate_clp

    !--------------------------------------------------------------------------
    ! Dimensional roll damping: L_p = C_l_p * (rho * V * S * b^2) / 4
    !--------------------------------------------------------------------------
    pure function dimensional_roll_damping(clp, rho, V, wing) result(Lp)
        real(dp), intent(in) :: clp, rho, V
        type(WingGeometry), intent(in) :: wing
        real(dp) :: Lp, S, b

        S = wing_area(wing)
        b = wing%span
        Lp = clp * (rho * V * S * b**2) / 4.0_dp
    end function dimensional_roll_damping

end module roll_damping_mod


module atmosphere_mod
    use constants
    implicit none

contains

    !--------------------------------------------------------------------------
    ! ISA Standard Atmosphere
    !--------------------------------------------------------------------------
    pure subroutine standard_atmosphere(altitude, T, p, rho)
        real(dp), intent(in) :: altitude
        real(dp), intent(out) :: T, p, rho

        real(dp), parameter :: T0 = 288.15_dp      ! Sea level temperature (K)
        real(dp), parameter :: p0 = 101325.0_dp    ! Sea level pressure (Pa)
        real(dp), parameter :: R = 287.0_dp        ! Gas constant (J/kg/K)
        real(dp), parameter :: g = 9.81_dp         ! Gravity (m/s^2)
        real(dp), parameter :: L = 0.0065_dp       ! Lapse rate (K/m)

        if (altitude < 11000.0_dp) then
            ! Troposphere
            T = T0 - L * altitude
            p = p0 * (T / T0)**5.2561_dp
        else
            ! Stratosphere (up to 20 km)
            T = 216.65_dp
            p = 22632.0_dp * exp(-g * (altitude - 11000.0_dp) / (R * T))
        end if

        rho = p / (R * T)
    end subroutine standard_atmosphere

end module atmosphere_mod


!==============================================================================
! Main Program
!==============================================================================
program roll_damping_calculator
    use constants
    use wing_geometry_mod
    use flow_regime_mod
    use lift_curve_slope_mod
    use roll_damping_mod
    use atmosphere_mod
    implicit none

    type(WingGeometry) :: wing_subsonic, wing_supersonic
    real(dp) :: mach, altitude, T, p, rho
    real(dp) :: clp, cl_alpha, mu, beta
    real(dp) :: clp_strip
    integer :: regime, i
    real(dp) :: mach_range(11)
    character(len=12) :: regime_str, le_type

    ! Mach number range for sweep
    mach_range = [0.3_dp, 0.6_dp, 0.8_dp, 0.9_dp, 1.0_dp, &
                  1.1_dp, 1.2_dp, 1.5_dp, 2.0_dp, 2.5_dp, 3.0_dp]

    write(*,'(A)') "======================================================================"
    write(*,'(A)') "Roll Damping Derivative (C_l_p) Calculator - FORTRAN"
    write(*,'(A)') "SUBSONIC and SUPERSONIC Regimes"
    write(*,'(A)') "======================================================================"

    !--------------------------------------------------------------------------
    ! Example 1: Subsonic - General Aviation Aircraft
    !--------------------------------------------------------------------------
    write(*,'(/A)') "======================================================================"
    write(*,'(A)')  "EXAMPLE 1: SUBSONIC - General Aviation Aircraft"
    write(*,'(A)')  "======================================================================"

    wing_subsonic = WingGeometry( &
        span = 11.0_dp, &
        root_chord = 1.63_dp, &
        tip_chord = 1.12_dp, &
        sweep_le = 2.0_dp * DEG2RAD &
    )

    mach = 0.2_dp
    altitude = 3000.0_dp
    call standard_atmosphere(altitude, T, p, rho)

    write(*,'(/A)') "Wing Geometry:"
    write(*,'(A,F8.2,A)') "  Wingspan:         ", wing_subsonic%span, " m"
    write(*,'(A,F8.3)')   "  Taper ratio:      ", taper_ratio(wing_subsonic)
    write(*,'(A,F8.2)')   "  Aspect ratio:     ", aspect_ratio(wing_subsonic)
    write(*,'(A,F8.1,A)') "  LE Sweep:         ", wing_subsonic%sweep_le * RAD2DEG, " deg"

    write(*,'(/A)') "Flight Conditions:"
    write(*,'(A,F8.2)')   "  Mach number:      ", mach
    write(*,'(A,F8.0,A)') "  Altitude:         ", altitude, " m"

    call calculate_clp(wing_subsonic, mach, clp, regime)
    cl_alpha = cl_alpha_3d(wing_subsonic, mach)

    write(*,'(/A,A,A)') "Results (", trim(regime_name(regime)), "):"
    write(*,'(A,F8.4,A)') "  C_L_alpha (3D):   ", cl_alpha, " 1/rad"
    write(*,'(A,F8.4,A)') "  C_l_p:            ", clp, " 1/rad"

    !--------------------------------------------------------------------------
    ! Example 2: Supersonic - Fighter Aircraft
    !--------------------------------------------------------------------------
    write(*,'(/A)') "======================================================================"
    write(*,'(A)')  "EXAMPLE 2: SUPERSONIC - Fighter Aircraft"
    write(*,'(A)')  "======================================================================"

    wing_supersonic = WingGeometry( &
        span = 9.45_dp, &
        root_chord = 5.0_dp, &
        tip_chord = 1.5_dp, &
        sweep_le = 40.0_dp * DEG2RAD &
    )

    mach = 1.6_dp
    altitude = 12000.0_dp
    call standard_atmosphere(altitude, T, p, rho)

    write(*,'(/A)') "Wing Geometry:"
    write(*,'(A,F8.2,A)') "  Wingspan:         ", wing_supersonic%span, " m"
    write(*,'(A,F8.3)')   "  Taper ratio:      ", taper_ratio(wing_supersonic)
    write(*,'(A,F8.2)')   "  Aspect ratio:     ", aspect_ratio(wing_supersonic)
    write(*,'(A,F8.1,A)') "  LE Sweep:         ", wing_supersonic%sweep_le * RAD2DEG, " deg"

    mu = mach_cone_angle(mach)
    if (is_subsonic_le(mach, wing_supersonic%sweep_le)) then
        le_type = "Subsonic LE"
    else
        le_type = "Supersonic LE"
    end if

    write(*,'(/A)') "Flight Conditions:"
    write(*,'(A,F8.2)')   "  Mach number:      ", mach
    write(*,'(A,F8.0,A)') "  Altitude:         ", altitude, " m"
    write(*,'(A,F8.1,A)') "  Mach angle (mu):  ", mu * RAD2DEG, " deg"
    write(*,'(A,A)')      "  Leading edge:     ", trim(le_type)

    call calculate_clp(wing_supersonic, mach, clp, regime)
    cl_alpha = cl_alpha_3d(wing_supersonic, mach)
    beta = compressibility_factor(mach)

    write(*,'(/A,A,A)') "Results (", trim(regime_name(regime)), "):"
    write(*,'(A,F8.4)')   "  beta = sqrt(M^2-1): ", beta
    write(*,'(A,F8.4,A)') "  C_L_alpha (3D):   ", cl_alpha, " 1/rad"
    write(*,'(A,F8.4,A)') "  C_l_p:            ", clp, " 1/rad"

    ! Strip theory comparison
    clp_strip = clp_supersonic_strip(wing_supersonic, mach, 100)
    write(*,'(A,F8.4,A)') "  C_l_p (strip):    ", clp_strip, " 1/rad"

    !--------------------------------------------------------------------------
    ! Example 3: Mach Sweep
    !--------------------------------------------------------------------------
    write(*,'(/A)') "======================================================================"
    write(*,'(A)')  "EXAMPLE 3: C_l_p vs MACH NUMBER"
    write(*,'(A)')  "======================================================================"

    write(*,'(/A,F5.2,A,F5.2,A,F4.0,A)') &
        "Wing: AR=", aspect_ratio(wing_supersonic), &
        ", lambda=", taper_ratio(wing_supersonic), &
        ", Sweep_LE=", wing_supersonic%sweep_le * RAD2DEG, " deg"

    write(*,'(A)') "--------------------------------------------------"
    write(*,'(A)') "  Mach |    Regime    |   C_L_a  |   C_l_p"
    write(*,'(A)') "--------------------------------------------------"

    do i = 1, 11
        mach = mach_range(i)
        call calculate_clp(wing_supersonic, mach, clp, regime)
        cl_alpha = cl_alpha_3d(wing_supersonic, mach)
        regime_str = regime_name(regime)

        write(*,'(F6.2,A,A12,A,F8.4,A,F8.4)') &
            mach, " | ", regime_str, " | ", cl_alpha, " | ", clp
    end do

    write(*,'(A)') "--------------------------------------------------"
    write(*,'(/A)') "Notes:"
    write(*,'(A)')  "  - Negative C_l_p indicates roll damping (stable)"
    write(*,'(A)')  "  - |C_l_p| generally decreases with increasing Mach"
    write(*,'(A)')  "  - Transonic region uses interpolation (simplified)"
    write(*,'(A)')  "  - Subsonic LE provides better supersonic damping"

end program roll_damping_calculator
