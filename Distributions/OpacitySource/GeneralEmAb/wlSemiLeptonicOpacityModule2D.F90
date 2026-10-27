!!! modified date: June 17, 2020
!!! Authors: Codes wirtten by Gang Guo, based on the formalism by
!!! Andreas Lohs & Gang Guo
!!! neutrino CC opacity is given by 2D integrals
!!! processes to be considered in this routine
!!! ReactionIndex=1: v + n -> e- + p  (opacity in 1/cm)
!!! ReactionIndex=2: vb + p -> e+ + n (opacity in 1/cm)
!!! ReactionIndex=3: vb + p + e- -> n (opacity in 1/cm)
!!! ReactionIndex=4: n -> vb + e- + p [emissivity in v/(s cm^3 MeV)]

!=======================================================================
!  wlSemiLeptonicOpacityModule2D   (front end)
!
!  Entry points for the semi-leptonic CC rates via 2D integrals. The
!  integration is done by one of three interchangeable methods in
!  IntegrationMethods/ (same physics, wlSemiLeptonicPhysicsModule2D):
!
!    iEmAb2D_Energy   = 1  ENERGY  "Standard": linear Gauss-Legendre in
!                          the nucleon energy E2 and lepton energy E3, the
!                          layout of the original GSI routines
!                          (nE x nE nodes).
!    iEmAb2D_Momentum = 2  MOMENTUM (DEFAULT): integration in the nucleon
!                          and lepton momenta with the domains split at
!                          thresholds/kinks (nE outer nodes per segment,
!                          MAX(6,MIN(12,nE/4)) inner nodes per segment).
!                          Fastest at a given accuracy; recommended.
!    iEmAb2D_Tail     = 3  TAIL: kink/Fermi-split segments plus a
!                          Gauss-Laguerre tail, no thermal cutoff; most
!                          accurate, ~3-5x the cost of MOMENTUM.
!
!  Calling sequence (the quadrature nodes are computed only once):
!    CALL InitOpacity_CC_2D( nE [, Method_Option] )   ! host, once
!    CALL Opacity_CC_2D( WhichCorrection, ReactionIndex, Enu, OpaA, &
!                        T, Mul, Mun, Mup, ml, mn, mp, Un, Up )
!                                                   ! host or device
!    CALL FinalizeOpacity_CC_2D()
!  InitOpacity_CC_2D (wlSemiLeptonicNodes2D) stores the method and its
!  Gauss-Legendre rules in module tables and copies them to the device
!  (OpenMP offload / OpenACC), so Opacity_CC_2D can be called inside
!  offloaded loops without further data movement. The default method is
!  iEmAb2D_Default (MOMENTUM unless built with EMAB2D_METHOD=ENERGY|TAIL);
!  Method_Option overrides it at run time. The meaning of nE depends on
!  the method (see above). Opacity_CC_2D returns 0 if Init was not called.
!
!  The stateless kernels Opacity_CC_2D_Energy_Nodes,
!  Opacity_CC_2D_Momentum_Nodes and Opacity_CC_2D_Tail_Nodes take the
!  quadrature rules as arguments instead, for callers that manage their
!  own nodes.
!=======================================================================
MODULE wlSemiLeptonicOpacityModule2D

  USE wlKindModule, ONLY: dp
  USE wlSemiLeptonicPhysicsModule2D, ONLY: gauleg, opacity_unit_factor
  USE wlSemiLeptonicIntegration2D_Energy, ONLY: &
    Opacity_CC_2D_Energy_Nodes
  USE wlSemiLeptonicIntegration2D_Momentum, ONLY: &
    Opacity_CC_2D_Momentum_Nodes, Inner_Rule
  USE wlSemiLeptonicIntegration2D_Tail, ONLY: &
    Opacity_CC_2D_Tail_Nodes
  USE wlSemiLeptonicNodes2D, ONLY: &
    InitOpacity_CC_2D, FinalizeOpacity_CC_2D, &
    iEmAb2D_Energy, iEmAb2D_Momentum, iEmAb2D_Tail, iEmAb2D_Default, &
    Initialized_CC2D, iMethod_CC2D, nE_CC2D, &
    nE2_CC2D, x2_CC2D, w2_CC2D, nE3_CC2D, x3_CC2D, w3_CC2D

  IMPLICIT NONE
  PRIVATE

  !===========================================================================
  !  Semi-leptonic charged-current (CC) neutrino rates via 2D integration
  !
  !  ReactionIndex:
  !    1 : nu capture on neutron
  !    2 : anti-nu capture on proton
  !    3 : inverse neutron decay
  !    4 : neutron decay emissivity
  !
  !  WhichCorrection:
  !    0 : LO only
  !    1 : LO + weak magnetism (WM)
  !    2 : LO + WM + pseudoscalar (PS)
  !    3 : LO + WM + PS + form-factor dependence
  !
  !  Notes:
  !   - For ReactionIndex=4, final-state neutrino blocking is not included.
  !===========================================================================

  !========================================================================
! Routine to calculate the neutrino CC opacity via 2D integrals
! Inputs:
! (1) WhichCorrection=0,1,2,3 for LO, LO+weak magnetism(WM),
! LO+WM+pseduoscalar(PS), LO+WM+PS+form factor dependencies
! (2) ReactionIndex=ReactionIndex=1,2,3,4 for nu capture on n, nub capture on p, inverse
! neutron decay, neutron decay emissivity
! (3) Enu: neutrino energy in MeV
! (4) xT: Temperature in MeV
! (5) xMul, xml: relativistic chemical potentials & masses of
! e^-(mu^-) in MeV
! (6) xMun,xMup,xmn,xmp,xUn,xUp: relativistic chemical
! potentials, masses, & potentials of neutron and proton in MeV
! relativsitic disperstion relation for nucleons,
! E2=SQRT(p^2+Mass2^2)+U2
! Outputs:
! (7) OpaA: opacity in cm^-1 (output_in_km_inverse = .false.); for neutron decay, emissivities
! in (s cm^3 MeV)^-1
!     **** blocking factor of final state neutrinos in neutron decay
!     is not considered ****
!=============================================================
  PUBLIC :: Opacity_CC_2D
  PUBLIC :: InitOpacity_CC_2D, FinalizeOpacity_CC_2D
  PUBLIC :: iEmAb2D_Energy, iEmAb2D_Momentum, iEmAb2D_Tail, iEmAb2D_Default
  PUBLIC :: Initialized_CC2D, iMethod_CC2D, nE_CC2D
  ! Stateless kernels (quadrature rules as arguments)
  PUBLIC :: Opacity_CC_2D_Energy_Nodes, Opacity_CC_2D_Momentum_Nodes, Opacity_CC_2D_Tail_Nodes
  PUBLIC :: Inner_Rule, gauleg, opacity_unit_factor

CONTAINS

SUBROUTINE Opacity_CC_2D(WhichCorrection, ReactionIndex, Enu, OpaA, xT, xMul, xMun, xMup, &
     xml, xmn, xmp, xUn, xUp)

  ! Method and quadrature rules from InitOpacity_CC_2D (module tables of
  ! wlSemiLeptonicNodes2D, also present on the device).
  INTEGER , INTENT(IN)  :: ReactionIndex, WhichCorrection
  REAL(DP), INTENT(IN)  :: Enu,xT,xMul,xMun,xMup,xml,xmn,xmp,xUn,xUp
  REAL(DP), INTENT(OUT) :: OpaA
#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET
#elif defined(WEAKLIB_OACC)
  !$ACC ROUTINE SEQ
#endif

  IF( .NOT. Initialized_CC2D ) THEN
    ! No I/O or STOP in device-callable code: InitOpacity_CC_2D is mandatory
    OpaA = 0.0d0
    RETURN
  END IF

  SELECT CASE( iMethod_CC2D )
  CASE( iEmAb2D_Energy )
    CALL Opacity_CC_2D_Energy_Nodes(WhichCorrection, ReactionIndex, Enu, OpaA, xT, xMul, xMun, xMup, &
         xml, xmn, xmp, xUn, xUp, nE2_CC2D, x2_CC2D, w2_CC2D, nE3_CC2D, x3_CC2D, w3_CC2D)
  CASE( iEmAb2D_Tail )
    CALL Opacity_CC_2D_Tail_Nodes(WhichCorrection, ReactionIndex, Enu, OpaA, xT, xMul, xMun, xMup, &
         xml, xmn, xmp, xUn, xUp, nE2_CC2D, x2_CC2D, w2_CC2D)
  CASE DEFAULT
    CALL Opacity_CC_2D_Momentum_Nodes(WhichCorrection, ReactionIndex, Enu, OpaA, xT, xMul, xMun, xMup, &
         xml, xmn, xmp, xUn, xUp, nE2_CC2D, x2_CC2D, w2_CC2D, nE3_CC2D, x3_CC2D, w3_CC2D)
  END SELECT

END SUBROUTINE Opacity_CC_2D

END MODULE wlSemiLeptonicOpacityModule2D
