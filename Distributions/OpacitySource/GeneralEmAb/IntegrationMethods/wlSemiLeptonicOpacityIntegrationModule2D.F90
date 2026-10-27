!=======================================================================
!  wlSemiLeptonicOpacityIntegrationModule2D
!
!  Two-phase, GPU-oriented driver layer for the MOMENTUM integration
!  method (wlSemiLeptonicIntegration2D_Momentum) of the semi-leptonic CC
!  2D integrals. All physics (amplitudes, kinematic bounds, occupations)
!  lives in wlSemiLeptonicPhysicsModule2D and is reused here — this
!  module holds no duplicated formulas.
!
!  Reactions (same indices as Opacity_CC_2D):
!    1: nu + n -> l- + p        2: nubar + p -> l+ + n
!    3: nubar + p + l- -> n     4: n -> nubar + l- + p (emissivity)
!
!  Structure (designed so the (iE2,iE3) quadrature can be exposed as a
!  collapsed GPU loop with a sum reduction):
!
!    phase 0  CC2D_PrepareCall : once per (zone, Enu, reaction).
!             Maps physical inputs (Mul, Mun, Mup, ...) to the internal
!             particle-2/3/4 convention, computes the p2 range
!             (Range_pn / Range_pn_D) and its kink-free outer segments
!             (P2_Segments, at most nSeg2Max), the per-call amplitude
!             coefficients (Calc_AmpCoeffs) and the constant prefactor.
!             The inner (p3) rule is the momentum rule precomputed by
!             InitOpacity_CC_2D (x3M_CC2D, wlSemiLeptonicNodes2D), as in
!             Opacity_CC_2D; CC2D_SetInnerRule (optional, after
!             PrepareCall) replaces it by any rule with 1 <= nE3 <= nE3Seg.
!    phase 1  CC2D_RowSetup    : once per outer Gauss-Legendre node iE2 =
!             1..nE2_CC2D (nodes precomputed by InitOpacity_CC_2D). A row
!             holds that node in
!             every outer segment: E2, weight, f2 occupation, the E3
!             window (Ebounds / Ebounds_D) and its kink-free inner
!             segments (E3_Segments). The (outer, inner) segment pairs
!             with a non-empty window are packed in summation order.
!    phase 2  CC2D_NodeKernel  : once per (row, node iE3),
!             iE3 = 1..nE3Max, where this module's nE3Max =
!             nSeg2Max*nSegMax*nE3Seg (packed segment pair (iE3-1)/nE3+1,
!             inner node MOD(iE3-1,nE3)+1 with nE3 = CD%nE3 <= nE3Seg;
!             nodes beyond the used pairs return 0, so the loop length
!             is fixed and is the maximum for any inner rule). NOTE:
!             this nE3Max is the row length and differs from the inner
!             nE3Max of wlSemiLeptonicIntegration2D_Momentum.
!             Amplitude x occupations x weights; returns the raw
!             integrand contribution. Accumulate over all nodes and
!             finish with CC2D_ApplyPrefactor (constant factors + sign).
!
!  A typical GPU driver (per timestep, OpenMP offload flavour), after
!  CALL InitOpacity_CC_2D( nE ) on the host:
!
!    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(CD,RD,res,tot)
!    DO iZone = 1, nZones; DO iE = 1, nEnu
!      CALL CC2D_PrepareCall(...,(iZone,iE)..., CD)
!      tot = 0
!      DO iE2 = 1, nE2_CC2D      ! or: hoist RowSetup into its own kernel
!        CALL CC2D_RowSetup(CD, iE2, RD)
!        !$OMP SIMD REDUCTION(+:tot) PRIVATE(res)
!        DO iE3 = 1, nE3Max
!          CALL CC2D_NodeKernel(CD, RD, iE3, res)
!          tot = tot + res
!        END DO
!      END DO
!      CALL CC2D_ApplyPrefactor(CD, tot)
!      Opacity(iE,iZone) = tot
!    END DO; END DO
!
!  The four per-reaction wrappers below run the same two phases serially
!  on the CPU with the rules of InitOpacity_CC_2D and are bit-identical to
!  Opacity_CC_2D with the MOMENTUM method; with any CC2D_SetInnerRule the
!  same loop reproduces Opacity_CC_2D_Momentum_Nodes.
!
!  NOTE on conventions: the wrappers take PHYSICAL arguments (the same
!  list as Opacity_CC_2D: Mul, Mun, Mup, ml, massn, massp, Un, Up); the
!  mapping to particle 2/3/4 (including the anti-lepton chemical
!  potential -Mul for reaction 2) is done internally, so callers can no
!  longer get the per-reaction argument shuffling wrong.
!=======================================================================
MODULE wlSemiLeptonicOpacityIntegrationModule2D

  USE wlKindModule, ONLY: dp
  USE wlEosConstantsModule, ONLY: pi, Gw_MeV, Vud, hbarc, hbar
  USE wlSemiLeptonicPhysicsModule2D, ONLY: &
    AmpCoeffsType, Calc_AmpCoeffs, Calc_Ampsq, &
    Range_pn, Range_pn_D, Ebounds, Ebounds_D, FD, gauleg, &
    opacity_unit_factor
  USE wlSemiLeptonicIntegration2D_Momentum, ONLY: &
    E2_Node, E3_Segments, E3_Node, P2_Segments, nSegMax, nE3Seg, nSeg2Max, &
    xE3Seg, wE3Seg
  USE wlSemiLeptonicNodes2D, ONLY: &
    nE2_CC2D, x2_CC2D, w2_CC2D, nE3M_CC2D, x3M_CC2D, w3M_CC2D

  IMPLICIT NONE
  PRIVATE

  ! Packed (outer, inner) segment pairs per row, and the per-row length
  ! of the CC2D_NodeKernel loop
  INTEGER, PARAMETER :: nPairMax = nSeg2Max*nSegMax
  INTEGER, PARAMETER :: nE3Max   = nPairMax*nE3Seg

  ! ---- phase-0 output: everything fixed for one (zone, Enu, reaction) call
  TYPE, PUBLIC :: CC2D_CallData
    INTEGER  :: ReactionIndex = 0
    INTEGER  :: IncludeCorrections(3) = 0
    REAL(DP) :: Enu = 0.0d0, T = 0.0d0
    REAL(DP) :: Mass2 = 0.0d0, Mass3 = 0.0d0, Mass4 = 0.0d0
    REAL(DP) :: U2 = 0.0d0, U4 = 0.0d0
    REAL(DP) :: Mu2 = 0.0d0, Mu3 = 0.0d0, Mu4 = 0.0d0
    REAL(DP) :: p2min = 0.0d0, p2max = 0.0d0, p4max = 0.0d0
    INTEGER  :: nSeg2 = 0                                         ! outer segments
    REAL(DP) :: p2a(nSeg2Max) = 0.0d0, p2b(nSeg2Max) = 0.0d0
    LOGICAL  :: HasPhaseSpace = .FALSE.
    TYPE(AmpCoeffsType) :: Coeffs
    INTEGER  :: nE3 = nE3Seg                                      ! inner rule on [0,1]
    REAL(DP) :: x3(nE3Seg) = xE3Seg, w3(nE3Seg) = wE3Seg          ! (per inner segment)
  END TYPE CC2D_CallData

  ! ---- phase-1 output: one outer (p2) quadrature row = one Gauss-Legendre
  !      node in each outer segment k = 1..nSeg2
  TYPE, PUBLIC :: CC2D_RowData
    REAL(DP) :: E2(nSeg2Max) = 0.0d0, wE2(nSeg2Max) = 0.0d0
    REAL(DP) :: f2(nSeg2Max) = 0.0d0, b2(nSeg2Max) = 0.0d0       ! b2 = 1 - f2
    REAL(DP) :: E3_min(nSeg2Max) = 0.0d0, E3_max(nSeg2Max) = 0.0d0
    INTEGER  :: nSeg(nSeg2Max) = 0                                ! E3 segments
    REAL(DP) :: p3a(nSegMax,nSeg2Max) = 0.0d0                     ! in p3
    REAL(DP) :: p3b(nSegMax,nSeg2Max) = 0.0d0
    INTEGER  :: nPair = 0                                         ! packed pairs:
    INTEGER  :: kPair(nPairMax) = 0, sPair(nPairMax) = 0          ! outer k, inner seg
    LOGICAL  :: HasWindow = .FALSE.
  END TYPE CC2D_RowData

  PUBLIC :: CC2D_PrepareCall
  PUBLIC :: CC2D_SetInnerRule
  PUBLIC :: CC2D_RowSetup
  PUBLIC :: CC2D_NodeKernel
  PUBLIC :: CC2D_ApplyPrefactor

  PUBLIC :: NuAbsorptionOnNeutrons
  PUBLIC :: NuBarAbsorptionOnProtons
  PUBLIC :: InverseNeutronDecay
  PUBLIC :: NeutronEmissivity

  PUBLIC :: gauleg     ! re-exported from the physics module
  PUBLIC :: nE2_CC2D   ! number of outer nodes (CC2D_RowSetup loop), set by InitOpacity_CC_2D
  PUBLIC :: nE3Max   ! per-row loop length for CC2D_NodeKernel (see header)

CONTAINS

!-----------------------------------------------------------------------
!  Phase 0: per-call setup
!-----------------------------------------------------------------------
SUBROUTINE CC2D_PrepareCall( ReactionIndex, Enu, T, Mul, Mun, Mup, ml, massn, massp, &
                             Un, Up, IncludeCorrections, CD )

  INTEGER , INTENT(IN)  :: ReactionIndex
  REAL(DP), INTENT(IN)  :: Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up
  INTEGER , INTENT(IN)  :: IncludeCorrections(3)
  TYPE(CC2D_CallData), INTENT(OUT) :: CD
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: p2min, p2max, anti

  CD%ReactionIndex      = ReactionIndex
  CD%IncludeCorrections = IncludeCorrections
  CD%Enu                = Enu
  CD%T                  = T
  CD%Mass3              = ml

  ! Inner rule precomputed by InitOpacity_CC_2D (as in Opacity_CC_2D);
  ! tabulated nE3Seg-point rule if Init was not called
  IF( nE3M_CC2D >= 1 ) THEN
    CD%nE3 = nE3M_CC2D
    CD%x3  = x3M_CC2D
    CD%w3  = w3M_CC2D
  ELSE
    CD%nE3 = nE3Seg
    CD%x3  = xE3Seg
    CD%w3  = wE3Seg
  END IF

  ! Particle-2/3/4 mapping, identical to Opacity_CC_2D's dispatch
  SELECT CASE(ReactionIndex)
  CASE (1)           ! nu + n -> l- + p
    CD%Mass2 = massn;  CD%Mass4 = massp
    CD%U2    = Un;     CD%U4    = Up
    CD%Mu2   = Mun;    CD%Mu3   = Mul;    CD%Mu4 = Mup
    anti     = 1.0d0
  CASE (2)           ! nubar + p -> l+ + n  (anti-lepton: Mu3 = -Mul)
    CD%Mass2 = massp;  CD%Mass4 = massn
    CD%U2    = Up;     CD%U4    = Un
    CD%Mu2   = Mup;    CD%Mu3   = -Mul;   CD%Mu4 = Mun
    anti     = -1.0d0
  CASE (3, 4)        ! inverse neutron decay / neutron decay
    CD%Mass2 = massp;  CD%Mass4 = massn
    CD%U2    = Up;     CD%U4    = Un
    CD%Mu2   = Mup;    CD%Mu3   = Mul;    CD%Mu4 = Mun
    anti     = -1.0d0
  CASE DEFAULT
    CD%HasPhaseSpace = .FALSE.
    RETURN
  END SELECT

  ! p2 integration range (phase space)
  IF( ReactionIndex <= 2 ) THEN
    CALL Range_pn( Enu, T, CD%Mass2, CD%Mass3, CD%Mass4, CD%U2, CD%U4, CD%Mu2, p2min, p2max )
    CD%p4max = 0.0d0
  ELSE
    CALL Range_pn_D( Enu, T, CD%Mass2, CD%Mass3, CD%Mass4, CD%U2, CD%U4, CD%Mu2, CD%Mu4, &
                     ReactionIndex, p2min, p2max, CD%p4max )
  END IF
  CD%HasPhaseSpace = ( p2min <= p2max )
  IF( CD%HasPhaseSpace ) THEN
    CD%p2min = p2min
    CD%p2max = p2max
  ELSE
    CD%p2min = 0.0d0
    CD%p2max = 0.0d0
  END IF

  ! Per-call amplitude coefficients
  CALL Calc_AmpCoeffs( CD%Mass2, CD%Mass4, CD%U2, CD%U4, anti, IncludeCorrections, CD%Coeffs )

  ! Kink-free outer segments (same call as Integral_2D / Integral_2D_D)
  CD%nSeg2 = 0
  CD%p2a   = 0.0d0
  CD%p2b   = 0.0d0
  IF( CD%HasPhaseSpace ) &
    CALL P2_Segments( MERGE( -1.0d0, 1.0d0, ReactionIndex <= 2 ), Enu, CD%Mass2, CD%Mass3, &
                      CD%Mass4, CD%U2, CD%U4, p2min, p2max, CD%nSeg2, CD%p2a, CD%p2b )

END SUBROUTINE CC2D_PrepareCall

!-----------------------------------------------------------------------
!  Optional, after CC2D_PrepareCall: replace the inner (p3) rule by
!  x3/w3 (nE3 nodes on [0,1], applied in each kink-free inner segment),
!  as in Opacity_CC_2D_Momentum_Nodes. nE3 outside [1,nE3Seg] keeps the
!  tabulated nE3Seg-point default. The NodeKernel loop length nE3Max is
!  unchanged (nodes beyond nPair*nE3 return 0).
!-----------------------------------------------------------------------
SUBROUTINE CC2D_SetInnerRule( CD, nE3, x3, w3 )

  TYPE(CC2D_CallData), INTENT(INOUT) :: CD
  INTEGER            , INTENT(IN)    :: nE3
  REAL(DP)           , INTENT(IN)    :: x3(nE3), w3(nE3)
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  IF( nE3 >= 1 .AND. nE3 <= nE3Seg ) THEN
    CD%nE3        = nE3
    CD%x3         = 0.0d0
    CD%w3         = 0.0d0
    CD%x3(1:nE3)  = x3(1:nE3)
    CD%w3(1:nE3)  = w3(1:nE3)
  ELSE
    CD%nE3 = nE3Seg
    CD%x3  = xE3Seg
    CD%w3  = wE3Seg
  END IF

END SUBROUTINE CC2D_SetInnerRule

!-----------------------------------------------------------------------
!  Apply the constant prefactor to the accumulated node sum, using the
!  same operation sequence as Integral_2D / Integral_2D_D so results are
!  bit-identical to Opacity_CC_2D.
!-----------------------------------------------------------------------
SUBROUTINE CC2D_ApplyPrefactor( CD, tot )

  TYPE(CC2D_CallData), INTENT(IN)    :: CD
  REAL(DP)           , INTENT(INOUT) :: tot
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  SELECT CASE(CD%ReactionIndex)
  CASE (1, 2, 3)
    tot = tot*(Gw_MeV*Vud)**2/16.0d0/(pi**5)/hbarc
    tot = tot*opacity_unit_factor
  CASE (4)
    tot = tot*(Gw_MeV*Vud)**2/32.0d0/(pi**7)*CD%Enu**2/hbarc**3/hbar
  CASE DEFAULT
    tot = 0.0d0
  END SELECT

END SUBROUTINE CC2D_ApplyPrefactor

!-----------------------------------------------------------------------
!  Phase 1: per-row setup for outer node iE2 (1..nE2_CC2D, Gauss-Legendre
!  rule on [0,1] of InitOpacity_CC_2D, mapped to p2 in each outer segment
!  [p2a(k), p2b(k)])
!-----------------------------------------------------------------------
SUBROUTINE CC2D_RowSetup( CD, iE2, RD )

  TYPE(CC2D_CallData), INTENT(IN)  :: CD
  INTEGER            , INTENT(IN)  :: iE2
  TYPE(CC2D_RowData) , INTENT(OUT) :: RD
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: E2, P2, xq0_min, xq0_max
  INTEGER  :: k, iSeg

  RD%HasWindow = .FALSE.
  RD%E2        = 0.0d0
  RD%wE2       = 0.0d0
  RD%f2        = 0.0d0
  RD%b2        = 0.0d0
  RD%E3_min    = 0.0d0
  RD%E3_max    = 0.0d0
  RD%nSeg      = 0
  RD%p3a       = 0.0d0
  RD%p3b       = 0.0d0
  RD%nPair     = 0
  RD%kPair     = 0
  RD%sPair     = 0
  IF( .NOT. CD%HasPhaseSpace ) RETURN

  DO k = 1, CD%nSeg2

    CALL E2_Node( CD%p2a(k), CD%p2b(k), x2_CC2D(iE2), w2_CC2D(iE2), CD%Mass2, CD%U2, E2, P2, RD%wE2(k) )

    IF( CD%ReactionIndex <= 2 ) THEN
      CALL Ebounds( E2, CD%Enu, CD%T, P2, CD%Mass2, CD%Mass3, CD%Mass4, CD%U2, CD%U4, &
                    xq0_min, xq0_max )
      RD%E3_min(k) = CD%Enu - xq0_max
      RD%E3_max(k) = CD%Enu - xq0_min
    ELSE
      CALL Ebounds_D( E2, CD%Enu, CD%T, P2, CD%Mass2, CD%Mass3, CD%Mass4, CD%U2, CD%U4, &
                      CD%Mu3, CD%p4max, CD%ReactionIndex, xq0_min, xq0_max )
      RD%E3_min(k) = xq0_min - CD%Enu
      RD%E3_max(k) = xq0_max - CD%Enu
    END IF

    RD%E2(k) = E2
    RD%f2(k) = FD( E2, CD%Mu2, CD%T )
    RD%b2(k) = FD( CD%Mu2, E2, CD%T )   ! 1 - f2 without cancellation
    IF( RD%E3_max(k) - RD%E3_min(k) <= 0.0d0 ) CYCLE

    CALL E3_Segments( MERGE( -1.0d0, 1.0d0, CD%ReactionIndex <= 2 ), CD%Enu, E2, P2, &
                      CD%Mass3, CD%Mass4, CD%U4, RD%E3_min(k), RD%E3_max(k), RD%nSeg(k), &
                      RD%p3a(:,k), RD%p3b(:,k) )

    ! Pack the (outer, inner) segment pairs in the summation order of
    ! Integral_2D / Integral_2D_D
    DO iSeg = 1, RD%nSeg(k)
      RD%nPair = RD%nPair + 1
      RD%kPair(RD%nPair) = k
      RD%sPair(RD%nPair) = iSeg
    END DO

  END DO

  RD%HasWindow = ( RD%nPair > 0 )

END SUBROUTINE CC2D_RowSetup

!-----------------------------------------------------------------------
!  Phase 2: per-node kernel, row node iE3 = 1..nE3Max (this module's
!  nE3Max: packed segment pair x nE3Seg; node iE3 is decoded with the
!  inner rule size CD%nE3 <= nE3Seg); returns the raw integrand
!  contribution (0 for nodes beyond the used pairs). Sum over all nodes
!  in order, then apply CC2D_ApplyPrefactor.
!-----------------------------------------------------------------------
SUBROUTINE CC2D_NodeKernel( CD, RD, iE3, res )

  TYPE(CC2D_CallData), INTENT(IN)  :: CD
  TYPE(CC2D_RowData) , INTENT(IN)  :: RD
  INTEGER            , INTENT(IN)  :: iE3
  REAL(DP)           , INTENT(OUT) :: res
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: E2, E3, wE3, E4, xamp
  INTEGER  :: iPair, k

  res = 0.0d0
  IF( .NOT. RD%HasWindow ) RETURN
  iPair = ( iE3 - 1 )/CD%nE3 + 1
  IF( iPair > RD%nPair ) RETURN

  ! Outer segment k; node of its inner segment sPair (E3_Node indexing)
  k = RD%kPair(iPair)
  CALL E3_Node( RD%p3a(:,k), RD%p3b(:,k), ( RD%sPair(iPair) - 1 )*CD%nE3 + MOD( iE3 - 1, CD%nE3 ) + 1, &
                CD%nE3, CD%x3, CD%w3, CD%Mass3, E3, wE3 )
  E2 = RD%E2(k)

  IF( CD%ReactionIndex <= 2 ) THEN
    E4 = CD%Enu + E2 - E3
    CALL Calc_Ampsq( CD%Enu, E2,  E3, CD%Mass2, CD%Mass3, CD%Mass4, CD%U2, CD%U4, &
                     CD%ReactionIndex, CD%IncludeCorrections, CD%Coeffs, xamp )
  ELSE
    E4 = CD%Enu + E3 + E2
    CALL Calc_Ampsq( CD%Enu, E2, -E3, CD%Mass2, CD%Mass3, CD%Mass4, CD%U2, CD%U4, &
                     CD%ReactionIndex, CD%IncludeCorrections, CD%Coeffs, xamp )
  END IF

  ! 1 - f_k is evaluated as FD(mu_k, E_k, T) (no cancellation)
  SELECT CASE(CD%ReactionIndex)
  CASE (1, 2)   ! captures
    res = RD%f2(k)*FD( CD%Mu3, E3, CD%T )*FD( CD%Mu4, E4, CD%T )*xamp*RD%wE2(k)*wE3/CD%Enu**2
  CASE (3)      ! inverse neutron decay
    res = RD%f2(k)*FD( E3, CD%Mu3, CD%T )*FD( CD%Mu4, E4, CD%T )*xamp*RD%wE2(k)*wE3/CD%Enu**2
  CASE (4)      ! neutron decay (*neutrino blocking not included*)
    res = RD%b2(k)*FD( CD%Mu3, E3, CD%T )*FD( E4, CD%Mu4, CD%T )*xamp*RD%wE2(k)*wE3/CD%Enu**2
  END SELECT

END SUBROUTINE CC2D_NodeKernel

!-----------------------------------------------------------------------
!  CPU reference wrappers (bit-identical to Opacity_CC_2D with the
!  MOMENTUM method; quadrature rules from InitOpacity_CC_2D)
!-----------------------------------------------------------------------
SUBROUTINE NuAbsorptionOnNeutrons( Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up, &
    IncludeCorrections, tot )

  REAL(DP), INTENT(IN)  :: Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up
  INTEGER , INTENT(IN)  :: IncludeCorrections(3)
  REAL(DP), INTENT(OUT) :: tot
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  CALL CC2D_Reference( 1, Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up, &
                       IncludeCorrections, tot )

END SUBROUTINE NuAbsorptionOnNeutrons

SUBROUTINE NuBarAbsorptionOnProtons( Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up, &
    IncludeCorrections, tot )

  REAL(DP), INTENT(IN)  :: Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up
  INTEGER , INTENT(IN)  :: IncludeCorrections(3)
  REAL(DP), INTENT(OUT) :: tot
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  CALL CC2D_Reference( 2, Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up, &
                       IncludeCorrections, tot )

END SUBROUTINE NuBarAbsorptionOnProtons

SUBROUTINE InverseNeutronDecay( Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up, &
    IncludeCorrections, tot )

  REAL(DP), INTENT(IN)  :: Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up
  INTEGER , INTENT(IN)  :: IncludeCorrections(3)
  REAL(DP), INTENT(OUT) :: tot
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  CALL CC2D_Reference( 3, Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up, &
                       IncludeCorrections, tot )

END SUBROUTINE InverseNeutronDecay

SUBROUTINE NeutronEmissivity( Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up, &
    IncludeCorrections, tot )

  REAL(DP), INTENT(IN)  :: Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up
  INTEGER , INTENT(IN)  :: IncludeCorrections(3)
  REAL(DP), INTENT(OUT) :: tot
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  CALL CC2D_Reference( 4, Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up, &
                       IncludeCorrections, tot )

END SUBROUTINE NeutronEmissivity

SUBROUTINE CC2D_Reference( ReactionIndex, Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up, &
    IncludeCorrections, tot )

  INTEGER , INTENT(IN)  :: ReactionIndex
  REAL(DP), INTENT(IN)  :: Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up
  INTEGER , INTENT(IN)  :: IncludeCorrections(3)
  REAL(DP), INTENT(OUT) :: tot
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  TYPE(CC2D_CallData) :: CD
  TYPE(CC2D_RowData)  :: RD
  REAL(DP) :: res
  INTEGER  :: iE2, iE3

  tot = 0.0d0
  CALL CC2D_PrepareCall( ReactionIndex, Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up, &
                         IncludeCorrections, CD )
  IF( .NOT. CD%HasPhaseSpace ) RETURN

  DO iE2 = 1, nE2_CC2D
    CALL CC2D_RowSetup( CD, iE2, RD )
    IF( .NOT. RD%HasWindow ) CYCLE
    DO iE3 = 1, nE3Max
      CALL CC2D_NodeKernel( CD, RD, iE3, res )
      tot = tot + res
    END DO
  END DO

  CALL CC2D_ApplyPrefactor( CD, tot )

END SUBROUTINE CC2D_Reference

END MODULE wlSemiLeptonicOpacityIntegrationModule2D
