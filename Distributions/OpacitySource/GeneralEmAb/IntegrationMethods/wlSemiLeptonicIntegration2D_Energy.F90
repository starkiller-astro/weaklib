!=======================================================================
!  wlSemiLeptonicIntegration2D_Energy     ("STANDARD" integration method)
!
!  Energy integration of the semi-leptonic CC 2D rates: the layout of the
!  original GSI routines (OriginalGSIroutines/Opacity_SemiLep_2D_GSI.f90),
!  i.e. the "standard" method. Not the default: wlSemiLeptonicOpacityModule2D
!  uses the momentum method unless built with EMAB2D_METHOD=ENERGY.
!
!  The outer integral runs over the nucleon energy E2 on [E2_min,E2_max]
!  (from Range_pn / Range_pn_D) and the inner integral over the lepton
!  energy E3 on [E3_min,E3_max] (from Ebounds / Ebounds_D), each with nE
!  Gauss-Legendre nodes mapped linearly from [0,1]. This is the layout the
!  code used before the switch to p2 (outer) and kink-free p3 segments
!  (inner, wlSemiLeptonicIntegration2D_Momentum); kinks of the angle-integrated |M|^2 inside the E3
!  window are not resolved, so convergence in nE is slower.
!
!  No physics is duplicated here: amplitudes, kinematic bounds, Fermi-Dirac
!  factors and nodes come from wlSemiLeptonicPhysicsModule2D. Same
!  units and particle mapping as Opacity_CC_2D. The kernel is stateless:
!  Opacity_CC_2D_Energy_Nodes takes the Gauss-Legendre rules on [0,1] for
!  E2 (nE2 nodes) and E3 (nE3 nodes) as arguments; Opacity_CC_2D passes
!  the rules precomputed by InitOpacity_CC_2D (nE2 = nE3 = nE).
!=======================================================================
MODULE wlSemiLeptonicIntegration2D_Energy

  USE wlKindModule, ONLY: dp
  USE wlEosConstantsModule, ONLY: &
   pi, Gw_MeV, Vud, hbarc, hbar
  USE wlSemiLeptonicPhysicsModule2D, ONLY: &
   AmpCoeffsType, Calc_AmpCoeffs, Calc_Ampsq, &
   Range_pn, Range_pn_D, Ebounds, Ebounds_D, FD, &
   nQuadMax, opacity_unit_factor

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: Opacity_CC_2D_Energy_Nodes

CONTAINS

SUBROUTINE Opacity_CC_2D_Energy_Nodes(WhichCorrection, ReactionIndex, Enu, OpaA, xT, xMul, xMun, xMup, &
     xml, xmn, xmp, xUn, xUp, nE2, x2, w2, nE3, x3, w3)

  ! Energy method with explicit quadrature rules on [0,1] (e.g. from
  ! gauleg): x2/w2 (nE2 nodes) for the outer E2 integral and x3/w3 (nE3
  ! nodes) for the inner E3 integral. No fixed-size buffers are involved,
  ! so nE2/nE3 are not clamped.
  INTEGER , INTENT(IN)  :: ReactionIndex, WhichCorrection, nE2, nE3
  REAL(DP), INTENT(IN)  :: Enu,xT,xMul,xMun,xMup,xml,xmn,xmp,xUn,xUp
  REAL(DP), INTENT(IN)  :: x2(nE2), w2(nE2), x3(nE3), w3(nE3)
  REAL(DP), INTENT(OUT) :: OpaA
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: anti
  INTEGER  :: IncludeCorrections(3)

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
     OpaA = 0.0d0
     RETURN
  END SELECT

  SELECT CASE(ReactionIndex)
  CASE (1)
    anti = 1.0d0
    CALL Integral_2D_Energy(Enu, xT, xmn, xml, xmp, xUn, xUp, xMun, xMul, xMup, anti, 1, IncludeCorrections, &
                            nE2, x2, w2, nE3, x3, w3, OpaA)
  CASE (2)
    anti = -1.0d0
    CALL Integral_2D_Energy(Enu, xT, xmp, xml, xmn, xUp, xUn, xMup, -xMul, xMun, anti, 2, IncludeCorrections, &
                            nE2, x2, w2, nE3, x3, w3, OpaA)
  CASE (3)
    anti = -1.0d0
    CALL Integral_2D_D_Energy(Enu, xT, xmp, xml, xmn, xUp, xUn, xMup, xMul, xMun, anti, 3, IncludeCorrections, &
                              nE2, x2, w2, nE3, x3, w3, OpaA)
  CASE (4)
    anti = -1.0d0
    CALL Integral_2D_D_Energy(Enu, xT, xmp, xml, xmn, xUp, xUn, xMup, xMul, xMun, anti, 4, IncludeCorrections, &
                              nE2, x2, w2, nE3, x3, w3, OpaA)
  CASE DEFAULT
    OpaA = 0.0d0
  END SELECT

END SUBROUTINE Opacity_CC_2D_Energy_Nodes

! Captures (ReactionIndex = 1, 2): linear nodes in E2 (outer) and E3 (inner)
SUBROUTINE Integral_2D_Energy(Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu3, Mu4, anti, ReactionIndex, IncludeCorrections, &
                              nE2, x2, w2, nE3, x3, w3, res)

  REAL(DP), INTENT(IN)  :: Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu3, Mu4, anti
  INTEGER , INTENT(IN)  :: ReactionIndex, IncludeCorrections(3), nE2, nE3
  REAL(DP), INTENT(IN)  :: x2(nE2), w2(nE2)   ! E2 Gauss-Legendre nodes/weights on [0,1]
  REAL(DP), INTENT(IN)  :: x3(nE3), w3(nE3)   ! E3 Gauss-Legendre nodes/weights on [0,1]
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
  REAL(DP) :: xf2,xb3,xb4   ! f2, 1-f3, 1-f4
  REAL(DP) :: p2min, p2max
  TYPE(AmpCoeffsType) :: Coeffs

  res = 0.0d0
  CALL Range_pn(Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, p2min, p2max)
  IF( p2min > p2max ) RETURN
  E2_min = SQRT(p2min**2 + Mass2**2) + U2
  E2_max = SQRT(p2max**2 + Mass2**2) + U2

  CALL Calc_AmpCoeffs( Mass2, Mass4, U2, U4, anti, IncludeCorrections, Coeffs )

  ! Loop over incoming particle (i.e. particle 2)
  DO i=1, nE2
    E2 = x2(i)*(E2_max - E2_min) + E2_min
    wE2 = w2(i)*(E2_max - E2_min)
    P2 = SQRT( MAX( (E2 - U2)**2 - Mass2**2, 0.0d0 ) )
    E2 = SQRT( P2**2 + Mass2**2 ) + U2   ! on-shell E2 consistent with P2
    CALL Ebounds( E2, Enu, T, P2, Mass2, Mass3, Mass4, U2, U4, xq0_min, xq0_max )
    E3_min = Enu - xq0_max
    E3_max = Enu - xq0_min
    IF( E3_max - E3_min <= 0.0d0 ) CYCLE   ! empty kinematic window at this E2
    xf2 = FD(E2, Mu2, T)
  ! Loop over outgoing particle (i.e. particle 3)
    DO j=1,nE3
      E3 = x3(j)*(E3_max - E3_min) + E3_min
      wE3 = w3(j)*(E3_max - E3_min)
      E4 = Enu + E2 - E3
      CALL Calc_Ampsq( Enu, E2, E3, Mass2, Mass3, Mass4, U2, U4, ReactionIndex, IncludeCorrections, Coeffs, xamp)
      xb3 = FD(Mu3, E3, T)
      xb4 = FD(Mu4, E4, T)
      res = res + xf2*xb3*xb4*xamp*wE2*wE3/Enu**2
    END DO
  END DO

  ! This is in 1/cm (hbarc in MeV cm)
  res = res*(Gw_MeV*Vud)**2/16.0d0/(pi**5)/hbarc
  ! opacity_unit_factor = 1 unless output_in_km_inverse is set
  res = res*opacity_unit_factor

END SUBROUTINE Integral_2D_Energy

! Inverse decay / decay (ReactionIndex = 3, 4): linear nodes in E2 and E3
SUBROUTINE Integral_2D_D_Energy(Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu3, Mu4, anti, ReactionIndex, IncludeCorrections, &
                              nE2, x2, w2, nE3, x3, w3, res)

  REAL(DP), INTENT(IN)  :: Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu3, Mu4, anti
  INTEGER , INTENT(IN)  :: ReactionIndex, IncludeCorrections(3), nE2, nE3
  REAL(DP), INTENT(IN)  :: x2(nE2), w2(nE2)   ! E2 Gauss-Legendre nodes/weights on [0,1]
  REAL(DP), INTENT(IN)  :: x3(nE3), w3(nE3)   ! E3 Gauss-Legendre nodes/weights on [0,1]
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
  REAL(DP) :: xf2,xf3,xf4,xb2,xb3,xb4   ! f_k and 1-f_k
  REAL(DP) :: p2min, p2max, p4max
  TYPE(AmpCoeffsType) :: Coeffs

  res = 0.0d0
  CALL Range_pn_D(Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu4, ReactionIndex, p2min, p2max, p4max)
  IF( p2min > p2max ) RETURN
  E2_min = SQRT(p2min**2 + Mass2**2) + U2
  E2_max = SQRT(p2max**2 + Mass2**2) + U2

  CALL Calc_AmpCoeffs( Mass2, Mass4, U2, U4, anti, IncludeCorrections, Coeffs )

  ! Loop over incoming particle (i.e. particle 2)
  DO i=1, nE2
    E2 = x2(i)*(E2_max - E2_min) + E2_min
    wE2 = w2(i)*(E2_max - E2_min)
    P2 = SQRT( MAX( (E2-U2)**2 - Mass2**2, 0.0d0 ) )
    E2 = SQRT( P2**2 + Mass2**2 ) + U2   ! on-shell E2 consistent with P2
    CALL Ebounds_D( E2, Enu, T, P2, Mass2, Mass3, Mass4, U2, U4, Mu3, p4max, ReactionIndex, xq0_min, xq0_max )
    E3_min = xq0_min - Enu
    E3_max = xq0_max - Enu
    IF( E3_max - E3_min <= 0.0d0 ) CYCLE   ! empty kinematic window at this E2
    xf2 = FD(E2, Mu2, T)
    xb2 = FD(Mu2, E2, T)
  ! Loop over outgoing particle (i.e. particle 3)
    DO j=1,nE3
      E3 = x3(j)*(E3_max - E3_min) + E3_min
      wE3 = w3(j)*(E3_max - E3_min)
      E4 = Enu + E3 + E2
      CALL Calc_Ampsq(Enu, E2, -E3, Mass2, Mass3, Mass4, U2, U4, ReactionIndex, IncludeCorrections, Coeffs, xamp)
      IF(ReactionIndex.eq.3) THEN ! inverse decay
        xf3 = FD(E3, Mu3, T)
        xb4 = FD(Mu4, E4, T)
        res = res + xf2*xf3*xb4*xamp*wE2*wE3/Enu**2
      ELSE IF(ReactionIndex.eq.4) THEN ! n decay, *blocking of neutrinos is not added*
        xb3 = FD(Mu3, E3, T)
        xf4 = FD(E4, Mu4, T)
        res = res + xb2*xb3*xf4*xamp*wE2*wE3/Enu**2
      END IF
    END DO
  END DO

  IF(ReactionIndex.eq.3) THEN
  ! This is in 1/cm (hbarc in MeV cm)
    res = res*(Gw_MeV*Vud)**2/16.0d0/(pi**5)/hbarc
  ! opacity_unit_factor = 1 unless output_in_km_inverse is set
    res = res*opacity_unit_factor
  ELSE IF(ReactionIndex.eq.4) THEN
    res = res*(Gw_MeV*Vud)**2/32.0d0/(pi**7)*Enu**2/hbarc**3/hbar
  END IF

END SUBROUTINE Integral_2D_D_Energy

END MODULE wlSemiLeptonicIntegration2D_Energy
