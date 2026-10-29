!!! modified date: June 17, 2020
!!! Authors: Codes wirtten by Gang Guo, based on the formalism by
!!! Andreas Lohs & Gang Guo
!!! neutrino CC opacity is given by 2D integrals
!!! processes to be considered in this routine
!!! ReactionIndex=1: v + n -> e- + p  (opacity in 1/cm)
!!! ReactionIndex=2: vb + p -> e+ + n (opacity in 1/cm)
!!! ReactionIndex=3: vb + p + e- -> n (opacity in 1/cm)
!!! ReactionIndex=4: n -> vb + e- + p [emissivity in v/(s cm^3 MeV)]

!!! nE defined in module: the No. of grids used in gaussian
!!! quadrature for each dimention
!!! seems nE needs be >~30 to reach an accuracy of 5%

MODULE wlSemiLeptonicOpacityModule2D

  USE wlKindModule, ONLY: dp
  USE wlEosConstantsModule, ONLY: pi, Gw_MeV, Vud, hbarc, hbar
  USE wlSemiLeptonicPhysicsModule2D, ONLY: &
    AmpCoeffsType, Calc_AmpCoeffs, Calc_Ampsq, &
    Range_pn, Range_pn_D, Ebounds, Ebounds_D, FD, FDPair, gauleg, &
    nQuadMax, opacity_unit_factor

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Opacity_CC_2D, Opacity_CC_2D_PrecomputedNodes
  PUBLIC :: AmpCoeffsType, Calc_AmpCoeffs, Calc_Ampsq
  PUBLIC :: Range_pn, Range_pn_D, Ebounds, Ebounds_D, FD, FDPair, gauleg
  PUBLIC :: opacity_unit_factor

CONTAINS

SUBROUTINE Opacity_CC_2D(WhichCorrection, ReactionIndex, Enu, OpaA, xT, xMul, xMun, xMup, &
     xml, xmn, xmp, xUn, xUp, nE)

  ! Back-compatible entry point: computes the [0,1] Gauss-Legendre nodes
  ! locally. In performance-critical (GPU / in-flight) code precompute the
  ! nodes once with gauleg and call Opacity_CC_2D_PrecomputedNodes instead.
  INTEGER , INTENT(IN)  :: ReactionIndex, WhichCorrection, nE
  REAL(DP), INTENT(IN)  :: Enu,xT,xMul,xMun,xMup,xml,xmn,xmp,xUn,xUp
  REAL(DP), INTENT(OUT) :: OpaA
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: xa(nQuadMax), wxa(nQuadMax)
  INTEGER  :: nEuse

  OpaA = 0.0d0
  IF( nE < 1 ) RETURN
  nEuse = MIN(nE, nQuadMax)
  CALL gauleg(0.d0, 1.d0, xa, wxa, nEuse)
  CALL Opacity_CC_2D_PrecomputedNodes(WhichCorrection, ReactionIndex, Enu, OpaA, xT, xMul, xMun, xMup, &
       xml, xmn, xmp, xUn, xUp, nEuse, xa, wxa)

END SUBROUTINE Opacity_CC_2D

SUBROUTINE Opacity_CC_2D_PrecomputedNodes(WhichCorrection, ReactionIndex, Enu, OpaA, xT, xMul, xMun, xMup, &
     xml, xmn, xmp, xUn, xUp, nE, xa, wxa)

  ! Same physics as Opacity_CC_2D, but takes precomputed Gauss-Legendre
  ! nodes/weights on [0,1] (from CALL gauleg(0.d0, 1.d0, xa, wxa, nE)).
  INTEGER , INTENT(IN)  :: ReactionIndex, WhichCorrection, nE
  REAL(DP), INTENT(IN)  :: Enu,xT,xMul,xMun,xMup,xml,xmn,xmp,xUn,xUp
  REAL(DP), INTENT(IN)  :: xa(nE), wxa(nE)
  REAL(DP), INTENT(OUT) :: OpaA
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: anti
  INTEGER  :: IncludeCorrections(3)

  OpaA = 0.0d0
  IF( nE < 1 .OR. nE > nQuadMax ) RETURN
  IF( Enu <= 0.0d0 .OR. xT <= 0.0d0 ) RETURN

  SELECT CASE(WhichCorrection)
  CASE (0)
     IncludeCorrections = (/ 0, 0, 0 /)
  CASE (1)
     IncludeCorrections = (/ 1, 0, 0 /)
  CASE (2)
     IncludeCorrections = (/ 1, 1, 0 /)
  CASE (3)
     IncludeCorrections = (/ 1, 1, 1 /)
  CASE DEFAULT
     ! Invalid input: no I/O or STOP in device-callable code; return zero
     RETURN
  END SELECT

  SELECT CASE(ReactionIndex)
  CASE (1)
    anti = 1.0d0
    CALL Integral_2D(Enu, xT, xmn, xml, xmp, xUn, xUp, xMun, xMul, xMup, anti, 1, IncludeCorrections, nE, xa, wxa, OpaA)
  CASE (2)
    anti = -1.0d0
    CALL Integral_2D(Enu, xT, xmp, xml, xmn, xUp, xUn, xMup, -xMul, xMun, anti, 2, IncludeCorrections, nE, xa, wxa, OpaA)
  CASE (3)
    anti = -1.0d0
    CALL Integral_2D_D(Enu, xT, xmp, xml, xmn, xUp, xUn, xMup, xMul, xMun, anti, 3, IncludeCorrections, nE, xa, wxa, OpaA)
  CASE (4)
    anti = -1.0d0
    CALL Integral_2D_D(Enu, xT, xmp, xml, xmn, xUp, xUn, xMup, xMul, xMun, anti, 4, IncludeCorrections, nE, xa, wxa, OpaA)
  CASE DEFAULT
    OpaA = 0.0d0
  END SELECT

END SUBROUTINE Opacity_CC_2D_PrecomputedNodes

SUBROUTINE Integral_2D &
  ( Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu3, Mu4, anti, &
    ReactionIndex, IncludeCorrections, nE, xa, wxa, res )

  REAL(DP), INTENT(IN)  :: Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu3, Mu4, anti
  INTEGER , INTENT(IN)  :: ReactionIndex, IncludeCorrections(3), nE
  REAL(DP), INTENT(IN)  :: xa(nE), wxa(nE)   ! Gauss-Legendre nodes/weights on [0,1]
  REAL(DP), INTENT(OUT) :: res
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: P2,xq0_min,xq0_max,E3_min,E3_max,E2_min,E2_max
  REAL(DP) :: wE2, wE3
  INTEGER  :: i, j
  REAL(DP) :: xamp, E2, E3, E4
  REAL(DP) :: xf2,xf3,xf4
  REAL(DP) :: p2min, p2max
  TYPE(AmpCoeffsType) :: Coeffs

  res = 0.0d0
  CALL Range_pn(Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, p2min, p2max)
  IF( p2min > p2max ) RETURN
  E2_min = SQRT(p2min**2 + Mass2**2) + U2
  E2_max = SQRT(p2max**2 + Mass2**2) + U2

  CALL Calc_AmpCoeffs( Mass2, Mass4, U2, U4, anti, IncludeCorrections, Coeffs )

  ! Loop over incoming particle (i.e. particle 2)
  DO i=1, nE
    E2 = xa(i)*(E2_max - E2_min) + E2_min
    wE2 = wxa(i)*(E2_max - E2_min)
    P2 = SQRT( MAX( 0.0d0, (E2 - U2)**2 - Mass2**2 ) )
    CALL Ebounds( E2, Enu, T, P2, Mass2, Mass3, Mass4, U2, U4, xq0_min, xq0_max )
    E3_min = Enu - xq0_max
    E3_max = Enu - xq0_min
    IF( E3_max - E3_min <= 0.0d0 ) CYCLE   ! empty kinematic window at this E2
    xf2 = FD(E2, Mu2, T)
  ! Loop over outgoing particle (i.e. particle 3)
    DO j=1,nE
      E3 = xa(j)*(E3_max - E3_min) + E3_min
      wE3 = wxa(j)*(E3_max - E3_min)
      E4 = Enu + E2 - E3
      CALL Calc_ampsq( Enu, E2, E3, Mass2, Mass3, Mass4, U2, U4, ReactionIndex, IncludeCorrections, Coeffs, xamp)
      IF( xamp < 0.0_DP )THEN
        res = -1.0_DP
        RETURN
      END IF
      ! FD(-E,-mu,T) is the vacancy, without subtracting from one.
      xf3 = FD(-E3, -Mu3, T)
      xf4 = FD(-E4, -Mu4, T)
      res = res + xf2*xf3*xf4*xamp*wE2*wE3/Enu**2
    END DO
  END DO

  ! This is in 1/cm
  res = res*(Gw_MeV*Vud)**2/16.0d0/(pi**5)/hbarc
  res = MAX( 0.0d0, res*opacity_unit_factor )

END SUBROUTINE Integral_2D

SUBROUTINE Integral_2D_D &
  ( Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu3, Mu4, anti, &
    ReactionIndex, IncludeCorrections, nE, xa, wxa, res )

  REAL(DP), INTENT(IN)  :: Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu3, Mu4, anti
  INTEGER , INTENT(IN)  :: ReactionIndex, IncludeCorrections(3), nE
  REAL(DP), INTENT(IN)  :: xa(nE), wxa(nE)   ! Gauss-Legendre nodes/weights on [0,1]
  REAL(DP), INTENT(OUT) :: res
#if defined(WEAKLIB_OMP_OL)
    !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
    !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: P2,xq0_min,xq0_max,E3_min,E3_max,E2_min,E2_max
  REAL(DP) :: wE2, wE3
  INTEGER :: i,j
  REAL(DP) :: xamp
  REAL(DP) :: E2, E3, E4
  REAL(DP) :: xf2,xf3,xf4, xb2,xb3,xb4
  REAL(DP) :: p2min, p2max, p4max
  TYPE(AmpCoeffsType) :: Coeffs

  res = 0.0d0
  CALL Range_pn_D(Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu4, ReactionIndex, p2min, p2max, p4max)
  IF( p2min > p2max ) RETURN
  E2_min = SQRT(p2min**2 + Mass2**2) + U2
  E2_max = SQRT(p2max**2 + Mass2**2) + U2

  CALL Calc_AmpCoeffs( Mass2, Mass4, U2, U4, anti, IncludeCorrections, Coeffs )

  ! Loop over incoming particle (i.e. particle 2)
  DO i=1, nE
    E2 = xa(i)*(E2_max - E2_min) + E2_min
    wE2 = wxa(i)*(E2_max - E2_min)
    P2 = SQRT( MAX( 0.0d0, (E2-U2)**2 - Mass2**2 ) )
    CALL Ebounds_D( E2, Enu, T, P2, Mass2, Mass3, Mass4, U2, U4, Mu3, p4max, ReactionIndex, xq0_min, xq0_max )
    E3_min = xq0_min - Enu
    E3_max = xq0_max - Enu
    IF( E3_max - E3_min <= 0.0d0 ) CYCLE   ! empty kinematic window at this E2
    CALL FDPair(E2, Mu2, T, xf2, xb2)
  ! Loop over outgoing particle (i.e. particle 3)
    DO j=1,nE
      E3 = xa(j)*(E3_max - E3_min) + E3_min
      wE3 = wxa(j)*(E3_max - E3_min)
      E4 = Enu + E3 + E2
      CALL Calc_ampsq(Enu, E2, -E3, Mass2, Mass3, Mass4, U2, U4, ReactionIndex, IncludeCorrections, Coeffs, xamp)
      IF( xamp < 0.0_DP )THEN
        res = -1.0_DP
        RETURN
      END IF
      CALL FDPair(E3, Mu3, T, xf3, xb3)
      CALL FDPair(E4, Mu4, T, xf4, xb4)
      IF(ReactionIndex.eq.3) THEN ! inverse decay
        res = res + xf2*xf3*xb4*xamp*wE2*wE3/Enu**2
      ELSE IF(ReactionIndex.eq.4) THEN ! n decay, *blocking of neutrinos is not added*
        res = res + xb2*xb3*xf4*xamp*wE2*wE3/Enu**2
      END IF
    END DO
  END DO

  IF(ReactionIndex.eq.3) THEN
  ! This is in 1/cm
    res = res*(Gw_MeV*Vud)**2/16.0d0/(pi**5)/hbarc
    res = MAX( 0.0d0, res*opacity_unit_factor )
  ELSE IF(ReactionIndex.eq.4) THEN
    res = res*(Gw_MeV*Vud)**2/32.0d0/(pi**7)*Enu**2/hbarc**3/hbar
  END IF

END SUBROUTINE Integral_2D_D

END MODULE wlSemiLeptonicOpacityModule2D
