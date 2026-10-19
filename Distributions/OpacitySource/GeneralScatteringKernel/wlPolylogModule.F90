!================================================================================
! Module: wlPolylogModule
!
! Implements the fast, analytically stabilized evaluation of Fermi-Dirac integrals
! F_0(z), F_1(z), F_2(z) using the polylogarithm expansions described in
! Chapter V, Section 8 of Bollig's Thesis (2018).
!================================================================================
MODULE wlPolylogModule

  USE wlKindModule, ONLY: DP

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Init_Polylogarithms, Fermi_Dirac_Integrals
  PUBLIC :: Set_Polylog_Backend

  ! Order of the Bernoulli expansion of S_2 and S_3 in powers of y = S_1(x).
  ! For |y| <= ln(2) (always the case with the z <= 0 / z > 0 split below),
  ! N_POLY = 18 gives a max. relative error of ~4e-16 for F_0, F_1, F_2 over
  ! z in [-50,50] (checked against mpmath); N_POLY = 8 gave ~7.5e-10.
  INTEGER, PARAMETER :: N_POLY = 18

  ! a_j = B_j/(j+1)! and c_j = C_3(j)/(j+1)! (Eqs. 8.9 - 8.12 of the thesis,
  ! B_1 = -1/2 convention), computed exactly with rational arithmetic.
  ! Hard-coded so that no initialization (and no mutable module state) is
  ! needed: thread safe and usable on the device.
  REAL(DP), PARAMETER :: a_poly(0:N_POLY) = [ &
     1.00000000000000000E+00_DP, -2.50000000000000000E-01_DP, &
     2.77777777777777762E-02_DP,  0.00000000000000000E+00_DP, &
    -2.77777777777777778E-04_DP,  0.00000000000000000E+00_DP, &
     4.72411186696900978E-06_DP,  0.00000000000000000E+00_DP, &
    -9.18577307466196408E-08_DP,  0.00000000000000000E+00_DP, &
     1.89788699889710005E-09_DP,  0.00000000000000000E+00_DP, &
    -4.06476164514422560E-11_DP,  0.00000000000000000E+00_DP, &
     8.92169102045645230E-13_DP,  0.00000000000000000E+00_DP, &
    -1.99392958607210744E-14_DP,  0.00000000000000000E+00_DP, &
     4.51898002961991825E-16_DP ]
  REAL(DP), PARAMETER :: c_poly(0:N_POLY) = [ &
     1.00000000000000000E+00_DP, -3.75000000000000000E-01_DP, &
     7.87037037037037063E-02_DP, -8.68055555555555594E-03_DP, &
     1.29629629629629631E-04_DP,  8.10185185185185162E-05_DP, &
    -3.41935716085375951E-06_DP, -1.32865646258503392E-06_DP, &
     8.66087175610985105E-08_DP,  2.52608759553203989E-08_DP, &
    -2.14469446836406487E-09_DP, -5.14011062201297940E-10_DP, &
     5.24958211460082972E-11_DP,  1.08877544066363182E-11_DP, &
    -1.27793960944936950E-12_DP, -2.36982417730874508E-13_DP, &
     3.10435788796546229E-14_DP,  5.26175862991250584E-15_DP, &
    -7.53847954994926526E-16_DP ]

CONTAINS

  !------------------------------------------------------------------------------
  ! Init_Polylogarithms
  ! No-op, kept for backward compatibility: the coefficients are parameters.
  !------------------------------------------------------------------------------
  SUBROUTINE Init_Polylogarithms()
  END SUBROUTINE Init_Polylogarithms


  ! This is in case in the future we want to implement different backends for the polylogarithm evaluation (e.g. direct numerical integration, or a different analytic expansion). For now, there is only one backend based on Bollig's thesis, so this is a no-op.
  SUBROUTINE Set_Polylog_Backend(iBackend)
    INTEGER, INTENT(IN) :: iBackend
  END SUBROUTINE Set_Polylog_Backend

  !------------------------------------------------------------------------------
  ! Fermi_Dirac_Integrals
  ! Evaluates F_0(z), F_1(z), F_2(z) analytically.
  ! Uses the mapping to Polylogarithms (Eq 8.3): F_n(z) = -n! S_{n+1}(-exp(z))
  ! and the inversion formulas (Eq 8.5) for z > 0, so that always
  ! |x| = exp(-|z|) <= 1 and |y| = |S_1(x)| <= ln(2).
  !------------------------------------------------------------------------------
  PURE SUBROUTINE Fermi_Dirac_Integrals( z, F0, F1, F2 )
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(IN)  :: z
    REAL(DP), INTENT(OUT) :: F0, F1, F2

    REAL(DP) :: x, y
    REAL(DP) :: S2, S3
    REAL(DP), PARAMETER :: pi_sq_over_6 = 1.644934066848226436_DP
    REAL(DP), PARAMETER :: pi_sq_over_3 = 3.289868133696452873_DP

    IF (z <= 0.0_DP) THEN
      ! --- nondegenerate: direct form ---
      x = -EXP(z)
      y = calc_S1(x)
      CALL eval_S2_S3(y, S2, S3)

      F0 = -y
      F1 = -S2
      F2 = -2.0_DP * S3
    ELSE
      ! --- degenerate: inversion form ---
      x = -EXP(-z)
      y = calc_S1(x)
      CALL eval_S2_S3(y, S2, S3)

      F0 = z - y
      F1 = 0.5_DP*z*z + pi_sq_over_6 + S2
      F2 = (1.0_DP/3.0_DP)*z**3 + pi_sq_over_3*z - 2.0_DP*S3
    END IF

  END SUBROUTINE Fermi_Dirac_Integrals


  !------------------------------------------------------------------------------
  ! calc_S1
  ! Implements Eq 8.6: Accelerated and precision-safe evaluation of S_1(x)
  ! S_1(x) = -ln(1-x)
  !------------------------------------------------------------------------------
  PURE FUNCTION calc_S1(x) RESULT(res)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(IN) :: x
    REAL(DP) :: res
    REAL(DP) :: abs_x

    abs_x = ABS(x)

    ! Protect against catastrophic cancellation for very small x
    IF (abs_x <= 4.0d-16) THEN
      res = x
    ELSE IF (abs_x <= 1.0d-4) THEN
      res = x + 0.5_DP*x**2 + (1.0_DP/3.0_DP)*x**3 + 0.25_DP*x**4
    ELSE
      res = -LOG(1.0_DP - x)
    END IF
  END FUNCTION calc_S1


  !------------------------------------------------------------------------------
  ! eval_S2_S3
  ! Evaluates the accelerated series expansions for S_2(x) and S_3(x) in terms
  ! of y = S_1(x) using Horner's method. (Eqs 8.9 and 8.10)
  !------------------------------------------------------------------------------
  PURE SUBROUTINE eval_S2_S3(y, S2, S3)
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif
    REAL(DP), INTENT(IN)  :: y
    REAL(DP), INTENT(OUT) :: S2, S3
    INTEGER :: j
    REAL(DP) :: term2, term3

    ! Start polynomial evaluation from highest order coefficient
    term2 = a_poly(N_POLY)
    term3 = c_poly(N_POLY)

    DO j = N_POLY - 1, 0, -1
      term2 = a_poly(j) + y * term2
      term3 = c_poly(j) + y * term3
    END DO

    ! The expansions are sum_{j=0} (coeff) * y^{j+1}
    S2 = y * term2
    S3 = y * term3
  END SUBROUTINE eval_S2_S3

END MODULE wlPolylogModule
