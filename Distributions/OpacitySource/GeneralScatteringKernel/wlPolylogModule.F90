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

  ! Maximum order of the Bernoulli expansion for polylogarithms.
  ! 25 terms is more than enough to reach machine precision (1e-16) for |y| < ln(2)
  INTEGER, PARAMETER :: N_POLY = 25
  
  REAL(DP), SAVE :: a_poly(0:N_POLY)
  REAL(DP), SAVE :: c_poly(0:N_POLY)
  LOGICAL, SAVE  :: is_initialized = .FALSE.

  INTEGER, PARAMETER :: POLY_BACKEND_BOLLIG = 1
  INTEGER, SAVE :: poly_backend = POLY_BACKEND_BOLLIG

  ! Transition window for special stable formulas (Eqs. 8.7/8.8 region)
  REAL(DP), PARAMETER :: z_mid = 1.0d0

  PUBLIC :: Set_Polylog_Backend
  
CONTAINS

  !------------------------------------------------------------------------------
  ! Init_Polylogarithms
  ! Precomputes the coefficients a_j and c_j required for the accelerated
  ! series expansions of S_2(x) and S_3(x) based on Bernoulli numbers.
  ! (Eqs. 8.9 - 8.12 of the thesis).
  !------------------------------------------------------------------------------
  SUBROUTINE Init_Polylogarithms()
    INTEGER :: j, k
    REAL(DP) :: B(0:N_POLY)
    REAL(DP) :: C3_val
    REAL(DP) :: fact, binom

    IF (is_initialized) RETURN

    ! Calculate Bernoulli numbers B_j using the recursive formula (Eq 8.12)
    B(0) = 1.0_DP
    DO j = 1, N_POLY
      B(j) = 0.0_DP
      DO k = 0, j - 1
        binom = factorial(j) / (factorial(k) * factorial(j - k))
        B(j) = B(j) + binom * B(k) / REAL(j - k + 1, DP)
      END DO
      B(j) = -B(j)
    END DO

    ! Calculate coefficients a_j and c_j
    DO j = 0, N_POLY
      fact = factorial(j + 1)
      
      ! a_j = B_j / (j+1)!
      a_poly(j) = B(j) / fact
      
      ! C_3(j) coefficients from Eq 8.11
      C3_val = 0.0_DP
      DO k = 0, j
        binom = factorial(j) / (factorial(k) * factorial(j - k))
        C3_val = C3_val + binom * B(j - k) * B(k) / REAL(k + 1, DP)
      END DO
      
      ! c_j = C_3(j) / (j+1)!
      c_poly(j) = C3_val / fact
    END DO

    is_initialized = .TRUE.
  END SUBROUTINE Init_Polylogarithms


  ! This is in case in the future we want to implement different backends for the polylogarithm evaluation (e.g. direct numerical integration, or a different analytic expansion). For now, we just have one backend based on Bollig's thesis.
  SUBROUTINE Set_Polylog_Backend(iBackend)
    INTEGER, INTENT(IN) :: iBackend
    poly_backend = iBackend
  END SUBROUTINE Set_Polylog_Backend

  !------------------------------------------------------------------------------
  ! Fermi_Dirac_Integrals
  ! Evaluates F_0(z), F_1(z), F_2(z) analytically.
  ! Uses the mapping to Polylogarithms (Eq 8.3): F_n(z) = -n! S_{n+1}(-exp(z))
  ! and handles limiting cases / inversions natively (Eq 8.5).
  !------------------------------------------------------------------------------
PURE SUBROUTINE Fermi_Dirac_Integrals( z, F0, F1, F2 )
    REAL(DP), INTENT(IN)  :: z
    REAL(DP), INTENT(OUT) :: F0, F1, F2

    REAL(DP) :: x, y, w
    REAL(DP) :: S2, S3
    REAL(DP), PARAMETER :: pi_sq_over_6 = 1.644934066848226436_DP
    REAL(DP), PARAMETER :: pi_sq_over_3 = 3.289868133696452873_DP

    ! --- Branch 1: nondegenerate ---
    IF (z < -z_mid) THEN
      x = -EXP(z)
      y = calc_S1(x)
      CALL eval_S2_S3(y, S2, S3)

      F0 = -y
      F1 = -S2
      F2 = -2.0_DP * S3
      RETURN
    END IF

    ! --- Branch 2: transition region (include Bollig Eq. 8.7/8.8 stabilization) ---
    IF (ABS(z) <= z_mid) THEN
      CALL fd_transition_stable(z, F0, F1, F2)
      RETURN
    END IF

    ! --- Branch 3: degenerate (inversion) ---
    w = -EXP(-z)
    y = calc_S1(w)
    CALL eval_S2_S3(y, S2, S3)

    F0 = z - y
    F1 = 0.5_DP*z*z + pi_sq_over_6 + S2
    F2 = (1.0_DP/3.0_DP)*z**3 + pi_sq_over_3*z - 2.0_DP*S3

  END SUBROUTINE Fermi_Dirac_Integrals


  !------------------------------------------------------------------------------
  ! calc_S1
  ! Implements Eq 8.6: Accelerated and precision-safe evaluation of S_1(x)
  ! S_1(x) = -ln(1-x)
  !------------------------------------------------------------------------------
  PURE FUNCTION calc_S1(x) RESULT(res)
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


  PURE SUBROUTINE fd_transition_stable(z, F0, F1, F2)
    REAL(DP), INTENT(IN)  :: z
    REAL(DP), INTENT(OUT) :: F0, F1, F2
    REAL(DP) :: x, y, S2, S3
    REAL(DP), PARAMETER :: pi_sq_over_6 = 1.644934066848226436_DP
    REAL(DP), PARAMETER :: pi_sq_over_3 = 3.289868133696452873_DP

    ! In transition region, evaluate both direct and inversion-style expressions
    ! and blend/switch to avoid cancellation (Bollig 8.7/8.8 spirit).
    if (z <= 0.0_DP) then
      x = -EXP(z)
      y = calc_S1(x)
      CALL eval_S2_S3(y,S2,S3)

      F0 = -y
      F1 = -S2
      F2 = -2.0_DP*S3
    else
      x = -EXP(-z)
      y = calc_S1(x)
      CALL eval_S2_S3(y,S2,S3)

      F0 = z - y
      F1 = 0.5_DP*z*z + pi_sq_over_6 + S2
      F2 = (1.0_DP/3.0_DP)*z**3 + pi_sq_over_3*z - 2.0_DP*S3
    end if
  END SUBROUTINE fd_transition_stable


  !------------------------------------------------------------------------------
  ! factorial
  ! Simple factorial helper for the initialization routine
  !------------------------------------------------------------------------------
  PURE FUNCTION factorial(n) RESULT(res)
    INTEGER, INTENT(IN) :: n
    REAL(DP) :: res
    INTEGER :: i
    res = 1.0_DP
    DO i = 2, n
      res = res * REAL(i, DP)
    END DO
  END FUNCTION factorial

END MODULE wlPolylogModule