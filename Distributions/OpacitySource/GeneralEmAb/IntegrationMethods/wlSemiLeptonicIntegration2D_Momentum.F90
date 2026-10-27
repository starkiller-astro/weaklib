!=======================================================================
!  wlSemiLeptonicIntegration2D_Momentum      (DEFAULT integration method)
!
!  Momentum integration of the semi-leptonic CC 2D rates: outer integral
!  in the nucleon momentum p2, split at the two-body thresholds of the
!  lepton window (P2_Segments); inner integral in the lepton momentum p3,
!  split at the kinks of the angle-integrated |M|^2 (E3_Segments). Both
!  variable changes remove SQRT endpoints of the integrand in energy, so
!  Gauss-Legendre converges much faster than in the energy method.
!  Physics: wlSemiLeptonicPhysicsModule2D. Thermal cutoff: Tfac.
!
!!! nE: the single resolution parameter. It is the No. of Gauss-Legendre
!!! nodes per segment of the outer (nucleon momentum p2) integral; the p2
!!! range is split at the two-body thresholds of the E3 window
!!! (P2_Segments, at most nSeg2Max segments, nE nodes each). The inner
!!! (lepton) integral is split at the kinks of the angle-integrated |M|^2
!!! and uses Inner_Count(nE) = MAX(nE3Min, MIN(nE3Seg, nE/4)) nodes per
!!! segment (in p3), see Inner_Rule. Opacity_CC_2D_Momentum_Nodes takes both rules
!!! explicitly (inner: 1..nE3Seg nodes).
!!! Accuracy (corrections 0-2, ThermoConditions.dat, 1-300 MeV, muons;
!!! error relative to the species' total opacity, species with total
!!! opacity < 1e-10 cm^-1 ignored), against a Tfac-free reference:
!!!   nE = 20: <~0.7% below E = max(mu_nu,0) + 8T,  nE = 24: <~0.7%,
!!!   nE = 32: <~0.03%,  nE = 40: <~0.01%,  nE >= 48: <~0.005%
!!! (strict pointwise max. rel. error at nE >= 48: ~5e-5, the Tfac floor).
!!! nE < 20 is not recommended. Correction 3 is limited to a few % above
!!! ~170 MeV by the discontinuous form-factor term, independently of nE.
!=======================================================================
MODULE wlSemiLeptonicIntegration2D_Momentum

  USE wlKindModule, ONLY: dp
  USE wlEosConstantsModule, ONLY: pi, Gw_MeV, Vud, hbarc, hbar
  USE wlSemiLeptonicPhysicsModule2D, ONLY: &
    AmpCoeffsType, Calc_AmpCoeffs, Calc_Ampsq, &
    Range_pn, Range_pn_D, Ebounds, Ebounds_D, Collinear_Roots, FD, gauleg, &
    nQuadMax, opacity_unit_factor

  IMPLICIT NONE
  PRIVATE

  ! Inner (lepton) quadrature. At fixed E2 the E3 window is split at the
  ! kinks of the angle-integrated |M|^2 (E3_Segments: at most 4 kinks, so
  ! at most nSegMax segments) and each segment is integrated in p3 with a
  ! Gauss-Legendre rule on [0,1] of at most nE3Seg nodes, passed down to
  ! Integral_2D(_D) / E3_Node as (nE3, x3, w3): Opacity_CC_2D passes
  ! Inner_Rule(nE), precomputed once by InitOpacity_CC_2D, and
  ! Opacity_CC_2D_Momentum_Nodes accepts any rule with 1 <= nE3 <= nE3Seg
  ! nodes. Loops over the inner nodes have the fixed length nE3Max; nodes
  ! of unused segments contribute nothing.
  INTEGER, PARAMETER :: nE3Seg  = 12
  INTEGER, PARAMETER :: nE3Min  = 6    ! fewest inner nodes used by Inner_Rule
  INTEGER, PARAMETER :: nSegMax = 5
  INTEGER, PARAMETER :: nE3Max  = nSegMax*nE3Seg
  REAL(DP), PARAMETER :: tGL12(6) = (/ &
    0.1252334085114689154724414d0, 0.3678314989981801937526915d0, &
    0.5873179542866174472967024d0, 0.7699026741943046870368938d0, &
    0.9041172563704748566784659d0, 0.9815606342467192506905491d0 /)
  REAL(DP), PARAMETER :: wGL12(6) = (/ &
    0.2491470458134027850005624d0, 0.2334925365383548087608499d0, &
    0.2031674267230659217490645d0, 0.1600783285433462263346525d0, &
    0.1069393259953184309602547d0, 0.0471753363865118271946160d0 /)
  REAL(DP), PARAMETER :: xE3Seg(nE3Seg) = (/ &
    0.5d0*(1.0d0-tGL12(6:1:-1)), 0.5d0*(1.0d0+tGL12(1:6)) /)
  REAL(DP), PARAMETER :: wE3Seg(nE3Seg) = (/ &
    0.5d0*wGL12(6:1:-1), 0.5d0*wGL12(1:6) /)

  ! Outer (p2) quadrature. [p2min, p2max] is split at the interior
  ! two-body thresholds of the E3 window (P2_Segments: at most nKinkMax
  ! points, so at most nSeg2Max segments) and each segment gets the nE
  ! caller nodes. Over the thermo/energy grid of the unit tests at most 2
  ! interior points occur; further points (never seen) would be dropped,
  ! which only costs accuracy.
  INTEGER, PARAMETER :: nKinkMax = 2
  INTEGER, PARAMETER :: nSeg2Max = nKinkMax + 1

  PUBLIC :: Opacity_CC_2D_Momentum_Nodes
  PUBLIC :: Inner_Rule
  ! Building blocks for the two-phase driver and the Tail method
  PUBLIC :: E2_Node, E3_Segments, E3_Node, P2_Segments
  PUBLIC :: nKinkMax, nSeg2Max
  PUBLIC :: nE3Seg, nE3Min, nSegMax, nE3Max, xE3Seg, wE3Seg

CONTAINS

!=======================================================================
!  Inner_Rule: inner (p3) Gauss-Legendre rule on [0,1] that goes with
!  the outer resolution nE: nE3 = MAX(nE3Min, MIN(nE3Seg, nE/4)) nodes
!  per kink-free segment (the tabulated xE3Seg/wE3Seg for nE3 = nE3Seg,
!  gauleg otherwise). The mapping was calibrated on muonic reactions:
!  fewer than 6 inner nodes fail at any nE, and the inner and outer
!  errors stay balanced along nE3 ~ nE/4.
!=======================================================================
SUBROUTINE Inner_Rule( nE, nE3, x3, w3 )

  INTEGER , INTENT(IN)  :: nE
  INTEGER , INTENT(OUT) :: nE3
  REAL(DP), INTENT(OUT) :: x3(nE3Seg), w3(nE3Seg)
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  nE3 = MAX( nE3Min, MIN( nE3Seg, nE/4 ) )
  IF( nE3 == nE3Seg ) THEN
    x3 = xE3Seg
    w3 = wE3Seg
  ELSE
    x3 = 0.0d0
    w3 = 0.0d0
    CALL gauleg( 0.d0, 1.d0, x3(1:nE3), w3(1:nE3), nE3 )
  END IF

END SUBROUTINE Inner_Rule

SUBROUTINE Opacity_CC_2D_Momentum_Nodes(WhichCorrection, ReactionIndex, Enu, OpaA, xT, xMul, xMun, xMup, &
     xml, xmn, xmp, xUn, xUp, nE2, x2, w2, nE3, x3, w3)

  ! Momentum method with explicit quadrature rules on [0,1] (e.g. from
  ! gauleg) for both integrals:
  !   x2/w2 (nE2): outer p2 rule, applied in each outer segment
  !   x3/w3 (nE3): inner p3 rule, applied in each kink-free inner
  !                segment, 1 <= nE3 <= nE3Seg (= 12)
  ! nE3 outside [1,nE3Seg] falls back to the tabulated nE3Seg-point rule
  ! (xE3Seg/wE3Seg). Opacity_CC_2D passes the rules precomputed by
  ! InitOpacity_CC_2D: gauleg(nE) outer, Inner_Rule(nE) inner.
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
  INTEGER  :: IncludeCorrections(3), nE3use
  REAL(DP) :: x3use(nE3Seg), w3use(nE3Seg)

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

  IF( nE3 >= 1 .AND. nE3 <= nE3Seg ) THEN
    nE3use = nE3
    x3use(1:nE3) = x3(1:nE3)
    w3use(1:nE3) = w3(1:nE3)
  ELSE
    nE3use = nE3Seg
    x3use  = xE3Seg
    w3use  = wE3Seg
  END IF

  SELECT CASE(ReactionIndex)
  CASE (1)
    anti = 1.0d0
    CALL Integral_2D(Enu, xT, xmn, xml, xmp, xUn, xUp, xMun, xMul, xMup, anti, 1, IncludeCorrections, &
                     nE2, x2, w2, nE3use, x3use, w3use, OpaA)
  CASE (2)
    anti = -1.0d0
    CALL Integral_2D(Enu, xT, xmp, xml, xmn, xUp, xUn, xMup, -xMul, xMun, anti, 2, IncludeCorrections, &
                     nE2, x2, w2, nE3use, x3use, w3use, OpaA)
  CASE (3)
    anti = -1.0d0
    CALL Integral_2D_D(Enu, xT, xmp, xml, xmn, xUp, xUn, xMup, xMul, xMun, anti, 3, IncludeCorrections, &
                       nE2, x2, w2, nE3use, x3use, w3use, OpaA)
  CASE (4)
    anti = -1.0d0
    CALL Integral_2D_D(Enu, xT, xmp, xml, xmn, xUp, xUn, xMup, xMul, xMun, anti, 4, IncludeCorrections, &
                       nE2, x2, w2, nE3use, x3use, w3use, OpaA)
  CASE DEFAULT
    OpaA = 0.0d0
  END SELECT

END SUBROUTINE Opacity_CC_2D_Momentum_Nodes

! For a more efficient parallelization this whole thing can be probably taken out, and we make 4 different integrations, one per reaction
SUBROUTINE Integral_2D(Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu3, Mu4, anti, ReactionIndex, IncludeCorrections, nE, xa, wxa, &
                       nE3, x3, w3, res)   ! captures

  REAL(DP), INTENT(IN)  :: Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu3, Mu4, anti
  INTEGER , INTENT(IN)  :: ReactionIndex, IncludeCorrections(3), nE, nE3
  REAL(DP), INTENT(IN)  :: xa(nE), wxa(nE)   ! outer (p2) Gauss-Legendre nodes/weights on [0,1]
  REAL(DP), INTENT(IN)  :: x3(nE3), w3(nE3)  ! inner (p3, per segment) rule on [0,1], nE3 <= nE3Seg
  REAL(DP), INTENT(OUT) :: res
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: P2,xq0_min,xq0_max,E3_min,E3_max
  REAL(DP) :: wE2, wE3
  INTEGER  :: i, j, nSeg
  REAL(DP) :: xamp, E2, E3, E4
  REAL(DP) :: xf2,xb3,xb4   ! f2, 1-f3, 1-f4
  REAL(DP) :: p2min, p2max
  REAL(DP) :: p3a(nSegMax), p3b(nSegMax)
  REAL(DP) :: p2a(nSeg2Max), p2b(nSeg2Max)
  INTEGER  :: k, nSeg2
  TYPE(AmpCoeffsType) :: Coeffs

  res = 0.0d0
  CALL Range_pn(Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, p2min, p2max)
  IF( p2min > p2max ) RETURN

  CALL Calc_AmpCoeffs( Mass2, Mass4, U2, U4, anti, IncludeCorrections, Coeffs )
  CALL P2_Segments( -1.0d0, Enu, Mass2, Mass3, Mass4, U2, U4, p2min, p2max, nSeg2, p2a, p2b )

  ! Loop over incoming particle (i.e. particle 2), integrated in p2:
  ! node i of each kink-free segment k (this order is reproduced by
  ! wlSemiLeptonicOpacityIntegrationModule2D)
  DO i=1, nE
  DO k=1, nSeg2
    CALL E2_Node( p2a(k), p2b(k), xa(i), wxa(i), Mass2, U2, E2, P2, wE2 )
    CALL Ebounds( E2, Enu, T, P2, Mass2, Mass3, Mass4, U2, U4, xq0_min, xq0_max )
    E3_min = Enu - xq0_max
    E3_max = Enu - xq0_min
    IF( E3_max - E3_min <= 0.0d0 ) CYCLE   ! empty kinematic window at this E2
    xf2 = FD(E2, Mu2, T)
    CALL E3_Segments( -1.0d0, Enu, E2, P2, Mass3, Mass4, U4, E3_min, E3_max, nSeg, p3a, p3b )
  ! Loop over outgoing particle (i.e. particle 3): kink-free segments in p3
    DO j=1, nSeg*nE3
      CALL E3_Node( p3a, p3b, j, nE3, x3, w3, Mass3, E3, wE3 )
      E4 = Enu + E2 - E3
      CALL Calc_ampsq( Enu, E2, E3, Mass2, Mass3, Mass4, U2, U4, ReactionIndex, IncludeCorrections, Coeffs, xamp)
      xb3 = FD(Mu3, E3, T)
      xb4 = FD(Mu4, E4, T)
      res = res + xf2*xb3*xb4*xamp*wE2*wE3/Enu**2
    END DO
  END DO
  END DO

  ! This is in 1/cm (hbarc in MeV cm)
  res = res*(Gw_MeV*Vud)**2/16.0d0/(pi**5)/hbarc
  ! opacity_unit_factor = 1 unless output_in_km_inverse is set
  res = res*opacity_unit_factor

END SUBROUTINE Integral_2D

! For a more efficient parallelization this whole thing can be probably taken out, and we make 4 different integrations, one per reaction
SUBROUTINE Integral_2D_D(Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu3, Mu4, anti, ReactionIndex, IncludeCorrections, nE, xa, wxa, &
                         nE3, x3, w3, res)   ! decay/inverse decay

  REAL(DP), INTENT(IN)  :: Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu3, Mu4, anti
  INTEGER , INTENT(IN)  :: ReactionIndex, IncludeCorrections(3), nE, nE3
  REAL(DP), INTENT(IN)  :: xa(nE), wxa(nE)   ! outer (p2) Gauss-Legendre nodes/weights on [0,1]
  REAL(DP), INTENT(IN)  :: x3(nE3), w3(nE3)  ! inner (p3, per segment) rule on [0,1], nE3 <= nE3Seg
  REAL(DP), INTENT(OUT) :: res
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: P2,xq0_min,xq0_max,E3_min,E3_max
  REAL(DP) :: wE2, wE3
  INTEGER :: i,j,nSeg
  REAL(DP) :: xamp
  REAL(DP) :: E2, E3, E4
  REAL(DP) :: xf2,xf3,xf4,xb2,xb3,xb4   ! f_k and 1-f_k
  REAL(DP) :: p2min, p2max, p4max
  REAL(DP) :: p3a(nSegMax), p3b(nSegMax)
  REAL(DP) :: p2a(nSeg2Max), p2b(nSeg2Max)
  INTEGER  :: k, nSeg2
  TYPE(AmpCoeffsType) :: Coeffs

  res = 0.0d0
  CALL Range_pn_D(Enu, T, Mass2, Mass3, Mass4, U2, U4, Mu2, Mu4, ReactionIndex, p2min, p2max, p4max)
  IF( p2min > p2max ) RETURN

  CALL Calc_AmpCoeffs( Mass2, Mass4, U2, U4, anti, IncludeCorrections, Coeffs )
  CALL P2_Segments( +1.0d0, Enu, Mass2, Mass3, Mass4, U2, U4, p2min, p2max, nSeg2, p2a, p2b )

  ! Loop over incoming particle (i.e. particle 2), integrated in p2:
  ! node i of each kink-free segment k (this order is reproduced by
  ! wlSemiLeptonicOpacityIntegrationModule2D)
  DO i=1, nE
  DO k=1, nSeg2
    CALL E2_Node( p2a(k), p2b(k), xa(i), wxa(i), Mass2, U2, E2, P2, wE2 )
    CALL Ebounds_D( E2, Enu, T, P2, Mass2, Mass3, Mass4, U2, U4, Mu3, p4max, ReactionIndex, xq0_min, xq0_max )
    E3_min = xq0_min - Enu
    E3_max = xq0_max - Enu
    IF( E3_max - E3_min <= 0.0d0 ) CYCLE   ! empty kinematic window at this E2
    xf2 = FD(E2, Mu2, T)
    xb2 = FD(Mu2, E2, T)
    CALL E3_Segments( +1.0d0, Enu, E2, P2, Mass3, Mass4, U4, E3_min, E3_max, nSeg, p3a, p3b )
  ! Loop over outgoing particle (i.e. particle 3): kink-free segments in p3
    DO j=1, nSeg*nE3
      CALL E3_Node( p3a, p3b, j, nE3, x3, w3, Mass3, E3, wE3 )
      E4 = Enu + E3 + E2
      CALL Calc_ampsq(Enu, E2, -E3, Mass2, Mass3, Mass4, U2, U4, ReactionIndex, IncludeCorrections, Coeffs, xamp)
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
  END DO

  IF(ReactionIndex.eq.3) THEN
  ! This is in 1/cm (hbarc in MeV cm)
    res = res*(Gw_MeV*Vud)**2/16.0d0/(pi**5)/hbarc
  ! opacity_unit_factor = 1 unless output_in_km_inverse is set
    res = res*opacity_unit_factor
  ELSE IF(ReactionIndex.eq.4) THEN
    res = res*(Gw_MeV*Vud)**2/32.0d0/(pi**7)*Enu**2/hbarc**3/hbar
  END IF

END SUBROUTINE Integral_2D_D

!=======================================================================
!  E2_Node: outer quadrature node, integrated in the nucleon momentum
!  p2 in [p2min, p2max] (dE2 = p2/E2* dp2). The angular integrals scale
!  with p2, so in E2 the integrand has a SQRT(E2 - E2min) endpoint that
!  limits Gauss-Legendre convergence; in p2 it is smooth.
!  x, w: Gauss-Legendre node and weight on [0,1].
!=======================================================================
SUBROUTINE E2_Node( p2min, p2max, x, w, Mass2, U2, E2, P2, wE2 )

  REAL(DP), INTENT(IN)  :: p2min, p2max, x, w, Mass2, U2
  REAL(DP), INTENT(OUT) :: E2, P2, wE2
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: E2f

  P2  = p2min + x*( p2max - p2min )
  E2f = SQRT( P2**2 + Mass2**2 )
  E2  = E2f + U2
  wE2 = w*( p2max - p2min )*P2/E2f

END SUBROUTINE E2_Node

!=======================================================================
!  P2_Segments: split the outer window [p2min, p2max] at the kinks of
!  the E3-integrated integrand. The dominant ones are the two-body
!  thresholds of the collinear configurations of Q0_Window_Envelope:
!  with W = Enu + E2 - U4 and P in { |Enu - P2|, Enu + P2 }, a pair of
!  collinear roots appears (Kallen lambda = 0) where
!    s(p2) = W**2 - P**2 = (Mass4 - s3*Mass3)**2
!  (captures, s3 = -1: s = (p3+p4)**2 >= (Mass4+Mass3)**2; decays,
!  s3 = +1: s = (p4-p3)**2 <= (Mass4-Mass3)**2; the other Kallen root is
!  unphysical). On the P = Enu + P2 branch this is the window threshold
!  (p2min of Range_pn / Range_pn_D); on the P = |Enu - P2| branch it lies
!  inside the window and the E3 window changes regime there: e.g. for
!  muonic inverse decay / decay at low Enu it sits next to the peak of
!  the p2 integrand and limited the outer convergence to ~n**-3.
!  With E2f = SQRT(P2**2+Mass2**2), a = Enu + U2 - U4, b = +/-1:
!    (a + E2f)**2 - (Enu + b*P2)**2 = mu**2
!    <=>  2*a*E2f = g + 2*b*Enu*P2,   g = mu**2 - a**2 - Mass2**2 + Enu**2,
!  a quadratic in P2 after squaring; the roots are returned
!  cancellation-free and the squared-in ones (a*(g+2*b*Enu*P2) <= 0)
!  rejected. Other regime changes (at-rest candidates, P2 = Enu, the
!  thermal q0 caps) were tested and gain less than their cost.
!  Returns nSeg2 >= 1 segments [p2a, p2b] that tile [p2min, p2max].
!=======================================================================
SUBROUTINE P2_Segments( s3, Enu, Mass2, Mass3, Mass4, U2, U4, p2min, p2max, nSeg2, p2a, p2b )

  REAL(DP), INTENT(IN)  :: s3, Enu, Mass2, Mass3, Mass4, U2, U4, p2min, p2max
  INTEGER , INTENT(OUT) :: nSeg2
  REAL(DP), INTENT(OUT) :: p2a(nSeg2Max), p2b(nSeg2Max)
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: a, g, qa, qb, qc, disc, q, r(2), b, eps, tmp, bp(nKinkMax+2)
  INTEGER  :: ib, ir, nb, i, j

  a   = Enu + U2 - U4
  g   = ( Mass4 - s3*Mass3 )**2 - a**2 - Mass2**2 + Enu**2
  qa  = ( a - Enu )*( a + Enu )
  qc  = ( a*Mass2 )**2 - 0.25d0*g**2
  eps = 1.0d-6*( p2max - p2min )

  ! Breakpoints: window edges plus interior thresholds, sorted
  nb = 1
  bp(1) = p2min
  DO ib = 1, 2
    b  = MERGE( 1.0d0, -1.0d0, ib == 1 )
    qb = -g*b*Enu
    ! qa*P2**2 + qb*P2 + qc = 0, discriminant a**2*(g**2 - 4*Mass2**2*qa)
    disc = a**2*( g**2 - 4.0d0*( Mass2**2 )*qa )
    IF( disc < 0.0d0 ) CYCLE
    q = -0.5d0*( qb + SIGN( SQRT(disc), qb ) )
    r = -1.0d0
    IF( qa /= 0.0d0 ) r(1) = q/qa
    IF( q  /= 0.0d0 ) r(2) = qc/q
    DO ir = 1, 2
      IF( a*( g + 2.0d0*b*Enu*r(ir) ) <= 0.0d0 ) CYCLE   ! squared-in root
      IF( r(ir) > p2min + eps .AND. r(ir) < p2max - eps .AND. nb < nKinkMax + 1 ) THEN
        nb = nb + 1
        bp(nb) = r(ir)
      END IF
    END DO
  END DO
  nb = nb + 1
  bp(nb) = p2max
  DO i = 3, nb - 1      ! insertion sort of the interior points
    tmp = bp(i)
    j   = i - 1
    DO WHILE( j >= 2 )
      IF( bp(j) <= tmp ) EXIT
      bp(j+1) = bp(j)
      j = j - 1
    END DO
    bp(j+1) = tmp
  END DO

  ! Drop (near-)duplicate points; segments tile [p2min, p2max] exactly
  nSeg2 = 0
  p2a   = 0.0d0
  p2b   = 0.0d0
  DO i = 1, nb - 1
    IF( bp(i+1) - bp(i) <= eps .AND. i < nb - 1 ) CYCLE
    nSeg2 = nSeg2 + 1
    p2a(nSeg2) = MERGE( p2min, p2b(MAX(nSeg2-1,1)), nSeg2 == 1 )
    p2b(nSeg2) = bp(i+1)
  END DO

END SUBROUTINE P2_Segments

!=======================================================================
!  E3_Segments: split the lepton window [E3_min, E3_max] at fixed E2 at
!  the kinks of the angle-integrated |M|^2. The kinks sit where the
!  min/max in the angular integration limits switch, i.e. at the
!  collinear configurations |p3 +/- p4| = P, P in {|Enu-P2|, Enu+P2}
!  (the same roots as in Q0_Window_Envelope; s3 = -1 capture, +1 decay).
!  Returns nSeg segments as p3 intervals [p3a, p3b] (integration in p3
!  removes the SQRT endpoint at E3 = Mass3).
!=======================================================================
SUBROUTINE E3_Segments( s3, Enu, E2, P2, Mass3, Mass4, U4, E3_min, E3_max, nSeg, p3a, p3b )

  REAL(DP), INTENT(IN)  :: s3, Enu, E2, P2, Mass3, Mass4, U4, E3_min, E3_max
  INTEGER , INTENT(OUT) :: nSeg
  REAL(DP), INTENT(OUT) :: p3a(nSegMax), p3b(nSegMax)
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: W, eps, tmp, y(2), bp(nSegMax+1)
  INTEGER  :: iP, ir, nRoot, nb, i, j

  W   = Enu - ( U4 - E2 )
  eps = 1.0d-9*( E3_max - E3_min )

  ! Breakpoints: window edges plus interior kinks, sorted
  nb = 1
  bp(1) = E3_min
  DO iP = 1, 2
    CALL Collinear_Roots( W, MERGE( ABS(Enu-P2), Enu+P2, iP == 1 ), Mass3, Mass4, nRoot, y )
    DO ir = 1, nRoot
      IF( y(ir) < Mass4 ) CYCLE
      tmp = s3*( y(ir) - W )
      IF( tmp > E3_min + eps .AND. tmp < E3_max - eps ) THEN
        nb = nb + 1
        bp(nb) = tmp
      END IF
    END DO
  END DO
  nb = nb + 1
  bp(nb) = E3_max
  DO i = 3, nb - 1      ! insertion sort of the interior points
    tmp = bp(i)
    j   = i - 1
    DO WHILE( j >= 2 )
      IF( bp(j) <= tmp ) EXIT
      bp(j+1) = bp(j)
      j = j - 1
    END DO
    bp(j+1) = tmp
  END DO

  nSeg = 0
  p3a  = 0.0d0
  p3b  = 0.0d0
  DO i = 1, nb - 1
    IF( bp(i+1) - bp(i) <= 1.0d-12*( ABS(E3_min) + ABS(E3_max) + 1.0d0 ) ) CYCLE
    nSeg = nSeg + 1
    p3a(nSeg) = SQRT( MAX( bp(i  )**2 - Mass3**2, 0.0d0 ) )
    p3b(nSeg) = SQRT( MAX( bp(i+1)**2 - Mass3**2, 0.0d0 ) )
  END DO

END SUBROUTINE E3_Segments

!=======================================================================
!  E3_Node: inner node j (1..nSeg*nE3) of the segments built by
!  E3_Segments, for the inner rule x3/w3 (nE3 nodes on [0,1]; the
!  default is nE3Seg, xE3Seg, wE3Seg): segment (j-1)/nE3+1, node
!  MOD(j-1,nE3)+1 in p3 (dE3 = p3/E3 dp3).
!=======================================================================
SUBROUTINE E3_Node( p3a, p3b, j, nE3, x3, w3, Mass3, E3, wE3 )

  REAL(DP), INTENT(IN)  :: p3a(nSegMax), p3b(nSegMax), Mass3
  INTEGER , INTENT(IN)  :: j, nE3
  REAL(DP), INTENT(IN)  :: x3(nE3), w3(nE3)
  REAL(DP), INTENT(OUT) :: E3, wE3
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: P3
  INTEGER  :: iSeg, iNode

  iSeg  = ( j - 1 )/nE3 + 1
  iNode = MOD( j - 1, nE3 ) + 1
  P3  = p3a(iSeg) + x3(iNode)*( p3b(iSeg) - p3a(iSeg) )
  E3  = SQRT( P3**2 + Mass3**2 )
  wE3 = w3(iNode)*( p3b(iSeg) - p3a(iSeg) )*P3/E3

END SUBROUTINE E3_Node

END MODULE wlSemiLeptonicIntegration2D_Momentum
