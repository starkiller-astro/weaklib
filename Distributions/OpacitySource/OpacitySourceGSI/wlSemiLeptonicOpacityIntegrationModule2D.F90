!=======================================================================
!  wlSemiLeptonicOpacityIntegrationModule2D
!
!  Two-phase, GPU-oriented driver layer for the semi-leptonic CC 2D
!  integrals. All physics (amplitudes, kinematic bounds, occupations)
!  lives in wlSemiLeptonicOpacityModule2D and is reused here — this
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
!             particle-2/3/4 convention, computes the E2 range
!             (Range_pn / Range_pn_D), the per-call amplitude
!             coefficients (Calc_AmpCoeffs) and the constant prefactor.
!    phase 1  CC2D_RowSetup    : once per E2 node.
!             E2, adjusted weight, f2 occupation and the E3 window
!             (Ebounds / Ebounds_D).
!    phase 2  CC2D_NodeKernel  : once per (E2,E3) node.
!             Amplitude x occupations x weights; returns the raw
!             integrand contribution. Accumulate over all nodes and
!             finish with CC2D_ApplyPrefactor (constant factors + sign).
!
!  A typical GPU driver (per timestep, OpenMP offload flavour):
!
!    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO COLLAPSE(2) PRIVATE(CD,RD,res,tot)
!    DO iZone = 1, nZones; DO iE = 1, nEnu
!      CALL CC2D_PrepareCall(...,(iZone,iE)..., CD)
!      tot = 0
!      DO iE2 = 1, nE            ! or: hoist RowSetup into its own kernel
!        CALL CC2D_RowSetup(CD, NodesE2(iE2), WeightsE2(iE2), RD)
!        !$OMP SIMD REDUCTION(+:tot) PRIVATE(res)
!        DO iE3 = 1, nE
!          CALL CC2D_NodeKernel(CD, RD, NodesE3(iE3), WeightsE3(iE3), res)
!          tot = tot + res
!        END DO
!      END DO
!      CALL CC2D_ApplyPrefactor(CD, tot)
!      Opacity(iE,iZone) = tot
!    END DO; END DO
!
!  The four per-reaction wrappers below run the same two phases serially
!  on the CPU and are bit-identical to Opacity_CC_2D(_PrecomputedNodes).
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
  USE wlSemiLeptonicOpacityModule2D, ONLY: &
    AmpCoeffsType, Calc_AmpCoeffs, Calc_Ampsq, &
    Range_pn, Range_pn_D, Ebounds, Ebounds_D, FD, gauleg, &
    opacity_unit_factor

  IMPLICIT NONE
  PRIVATE

  ! ---- phase-0 output: everything fixed for one (zone, Enu, reaction) call
  TYPE, PUBLIC :: CC2D_CallData
    INTEGER  :: ReactionIndex = 0
    INTEGER  :: IncludeCorrections(3) = 0
    REAL(DP) :: Enu = 0.0d0, T = 0.0d0
    REAL(DP) :: Mass2 = 0.0d0, Mass3 = 0.0d0, Mass4 = 0.0d0
    REAL(DP) :: U2 = 0.0d0, U4 = 0.0d0
    REAL(DP) :: Mu2 = 0.0d0, Mu3 = 0.0d0, Mu4 = 0.0d0
    REAL(DP) :: E2_min = 0.0d0, E2_max = 0.0d0, p4max = 0.0d0
    LOGICAL  :: HasPhaseSpace = .FALSE.
    TYPE(AmpCoeffsType) :: Coeffs
  END TYPE CC2D_CallData

  ! ---- phase-1 output: one E2 quadrature row
  TYPE, PUBLIC :: CC2D_RowData
    REAL(DP) :: E2 = 0.0d0, wE2 = 0.0d0, f2 = 0.0d0
    REAL(DP) :: E3_min = 0.0d0, E3_max = 0.0d0
    LOGICAL  :: HasWindow = .FALSE.
  END TYPE CC2D_RowData

  PUBLIC :: CC2D_PrepareCall
  PUBLIC :: CC2D_RowSetup
  PUBLIC :: CC2D_NodeKernel
  PUBLIC :: CC2D_ApplyPrefactor

  PUBLIC :: NuAbsorptionOnNeutrons
  PUBLIC :: NuBarAbsorptionOnProtons
  PUBLIC :: InverseNeutronDecay
  PUBLIC :: NeutronEmissivity

  PUBLIC :: gauleg   ! re-exported from wlSemiLeptonicOpacityModule2D

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

  ! E2 integration range (phase space)
  IF( ReactionIndex <= 2 ) THEN
    CALL Range_pn( Enu, T, CD%Mass2, CD%Mass3, CD%Mass4, CD%U2, CD%U4, CD%Mu2, p2min, p2max )
    CD%p4max = 0.0d0
  ELSE
    CALL Range_pn_D( Enu, T, CD%Mass2, CD%Mass3, CD%Mass4, CD%U2, CD%U4, CD%Mu2, CD%Mu4, &
                     ReactionIndex, p2min, p2max, CD%p4max )
  END IF
  CD%HasPhaseSpace = ( p2min <= p2max )
  IF( CD%HasPhaseSpace ) THEN
    CD%E2_min = SQRT(p2min**2 + CD%Mass2**2) + CD%U2
    CD%E2_max = SQRT(p2max**2 + CD%Mass2**2) + CD%U2
  ELSE
    CD%E2_min = 0.0d0
    CD%E2_max = 0.0d0
  END IF

  ! Per-call amplitude coefficients
  CALL Calc_AmpCoeffs( CD%Mass2, CD%Mass4, CD%U2, CD%U4, anti, IncludeCorrections, CD%Coeffs )

END SUBROUTINE CC2D_PrepareCall

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
!  Phase 1: per-E2-row setup (NodeE2/WeightE2: Gauss-Legendre on [0,1])
!-----------------------------------------------------------------------
SUBROUTINE CC2D_RowSetup( CD, NodeE2, WeightE2, RD )

  TYPE(CC2D_CallData), INTENT(IN)  :: CD
  REAL(DP)           , INTENT(IN)  :: NodeE2, WeightE2
  TYPE(CC2D_RowData) , INTENT(OUT) :: RD
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: E2, P2, xq0_min, xq0_max

  RD%HasWindow = .FALSE.
  RD%E3_min    = 0.0d0
  RD%E3_max    = 0.0d0
  IF( .NOT. CD%HasPhaseSpace ) RETURN

  E2     = NodeE2*(CD%E2_max - CD%E2_min) + CD%E2_min
  RD%wE2 = WeightE2*(CD%E2_max - CD%E2_min)
  P2     = SQRT( (E2 - CD%U2)**2 - CD%Mass2**2 )

  IF( CD%ReactionIndex <= 2 ) THEN
    CALL Ebounds( E2, CD%Enu, CD%T, P2, CD%Mass2, CD%Mass3, CD%Mass4, CD%U2, CD%U4, &
                  xq0_min, xq0_max )
    RD%E3_min = CD%Enu - xq0_max
    RD%E3_max = CD%Enu - xq0_min
  ELSE
    CALL Ebounds_D( E2, CD%Enu, CD%T, P2, CD%Mass2, CD%Mass3, CD%Mass4, CD%U2, CD%U4, &
                    CD%Mu3, CD%p4max, CD%ReactionIndex, xq0_min, xq0_max )
    RD%E3_min = xq0_min - CD%Enu
    RD%E3_max = xq0_max - CD%Enu
  END IF

  RD%E2        = E2
  RD%f2        = FD( E2, CD%Mu2, CD%T )
  RD%HasWindow = ( RD%E3_max - RD%E3_min > 0.0d0 )

END SUBROUTINE CC2D_RowSetup

!-----------------------------------------------------------------------
!  Phase 2: per-(E2,E3)-node kernel; returns the raw integrand
!  contribution. Sum over all nodes, then multiply by CD%Prefactor.
!-----------------------------------------------------------------------
SUBROUTINE CC2D_NodeKernel( CD, RD, NodeE3, WeightE3, res )

  TYPE(CC2D_CallData), INTENT(IN)  :: CD
  TYPE(CC2D_RowData) , INTENT(IN)  :: RD
  REAL(DP)           , INTENT(IN)  :: NodeE3, WeightE3
  REAL(DP)           , INTENT(OUT) :: res
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  REAL(DP) :: E3, wE3, E4, xamp, xf3, xf4

  res = 0.0d0
  IF( .NOT. RD%HasWindow ) RETURN

  E3  = NodeE3*(RD%E3_max - RD%E3_min) + RD%E3_min
  wE3 = WeightE3*(RD%E3_max - RD%E3_min)

  IF( CD%ReactionIndex <= 2 ) THEN
    E4 = CD%Enu + RD%E2 - E3
    CALL Calc_Ampsq( CD%Enu, RD%E2,  E3, CD%Mass2, CD%Mass3, CD%Mass4, CD%U2, CD%U4, &
                     CD%ReactionIndex, CD%IncludeCorrections, CD%Coeffs, xamp )
  ELSE
    E4 = CD%Enu + E3 + RD%E2
    CALL Calc_Ampsq( CD%Enu, RD%E2, -E3, CD%Mass2, CD%Mass3, CD%Mass4, CD%U2, CD%U4, &
                     CD%ReactionIndex, CD%IncludeCorrections, CD%Coeffs, xamp )
  END IF

  xf3 = FD( E3, CD%Mu3, CD%T )
  xf4 = FD( E4, CD%Mu4, CD%T )

  SELECT CASE(CD%ReactionIndex)
  CASE (1, 2)   ! captures
    res = RD%f2*(1.0d0-xf3)*(1.0d0-xf4)*xamp*RD%wE2*wE3/CD%Enu**2
  CASE (3)      ! inverse neutron decay
    res = RD%f2*xf3*(1.0d0-xf4)*xamp*RD%wE2*wE3/CD%Enu**2
  CASE (4)      ! neutron decay (*neutrino blocking not included*)
    res = (1.0d0-RD%f2)*(1.0d0-xf3)*xf4*xamp*RD%wE2*wE3/CD%Enu**2
  END SELECT

END SUBROUTINE CC2D_NodeKernel

!-----------------------------------------------------------------------
!  CPU reference wrappers (bit-identical to Opacity_CC_2D)
!-----------------------------------------------------------------------
SUBROUTINE NuAbsorptionOnNeutrons( Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up, &
    IncludeCorrections, nE, NodesE2, WeightsE2, NodesE3, WeightsE3, tot )

  REAL(DP), INTENT(IN)  :: Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up
  INTEGER , INTENT(IN)  :: IncludeCorrections(3), nE
  REAL(DP), INTENT(IN)  :: NodesE2(nE), WeightsE2(nE), NodesE3(nE), WeightsE3(nE)
  REAL(DP), INTENT(OUT) :: tot
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  CALL CC2D_Reference( 1, Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up, &
                       IncludeCorrections, nE, NodesE2, WeightsE2, NodesE3, WeightsE3, tot )

END SUBROUTINE NuAbsorptionOnNeutrons

SUBROUTINE NuBarAbsorptionOnProtons( Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up, &
    IncludeCorrections, nE, NodesE2, WeightsE2, NodesE3, WeightsE3, tot )

  REAL(DP), INTENT(IN)  :: Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up
  INTEGER , INTENT(IN)  :: IncludeCorrections(3), nE
  REAL(DP), INTENT(IN)  :: NodesE2(nE), WeightsE2(nE), NodesE3(nE), WeightsE3(nE)
  REAL(DP), INTENT(OUT) :: tot
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  CALL CC2D_Reference( 2, Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up, &
                       IncludeCorrections, nE, NodesE2, WeightsE2, NodesE3, WeightsE3, tot )

END SUBROUTINE NuBarAbsorptionOnProtons

SUBROUTINE InverseNeutronDecay( Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up, &
    IncludeCorrections, nE, NodesE2, WeightsE2, NodesE3, WeightsE3, tot )

  REAL(DP), INTENT(IN)  :: Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up
  INTEGER , INTENT(IN)  :: IncludeCorrections(3), nE
  REAL(DP), INTENT(IN)  :: NodesE2(nE), WeightsE2(nE), NodesE3(nE), WeightsE3(nE)
  REAL(DP), INTENT(OUT) :: tot
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  CALL CC2D_Reference( 3, Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up, &
                       IncludeCorrections, nE, NodesE2, WeightsE2, NodesE3, WeightsE3, tot )

END SUBROUTINE InverseNeutronDecay

SUBROUTINE NeutronEmissivity( Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up, &
    IncludeCorrections, nE, NodesE2, WeightsE2, NodesE3, WeightsE3, tot )

  REAL(DP), INTENT(IN)  :: Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up
  INTEGER , INTENT(IN)  :: IncludeCorrections(3), nE
  REAL(DP), INTENT(IN)  :: NodesE2(nE), WeightsE2(nE), NodesE3(nE), WeightsE3(nE)
  REAL(DP), INTENT(OUT) :: tot
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  CALL CC2D_Reference( 4, Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up, &
                       IncludeCorrections, nE, NodesE2, WeightsE2, NodesE3, WeightsE3, tot )

END SUBROUTINE NeutronEmissivity

SUBROUTINE CC2D_Reference( ReactionIndex, Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up, &
    IncludeCorrections, nE, NodesE2, WeightsE2, NodesE3, WeightsE3, tot )

  INTEGER , INTENT(IN)  :: ReactionIndex
  REAL(DP), INTENT(IN)  :: Enu, T, Mul, Mun, Mup, ml, massn, massp, Un, Up
  INTEGER , INTENT(IN)  :: IncludeCorrections(3), nE
  REAL(DP), INTENT(IN)  :: NodesE2(nE), WeightsE2(nE), NodesE3(nE), WeightsE3(nE)
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

  DO iE2 = 1, nE
    CALL CC2D_RowSetup( CD, NodesE2(iE2), WeightsE2(iE2), RD )
    IF( .NOT. RD%HasWindow ) CYCLE
    DO iE3 = 1, nE
      CALL CC2D_NodeKernel( CD, RD, NodesE3(iE3), WeightsE3(iE3), res )
      tot = tot + res
    END DO
  END DO

  CALL CC2D_ApplyPrefactor( CD, tot )

END SUBROUTINE CC2D_Reference

END MODULE wlSemiLeptonicOpacityIntegrationModule2D
