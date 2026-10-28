!=======================================================================
!  wlSemiLeptonicNodes2D
!
!  Quadrature state of the semi-leptonic CC 2D rates: the integration
!  method and its Gauss-Legendre rules on [0,1], computed ONCE by
!  InitOpacity_CC_2D on the host and copied to the device. Opacity_CC_2D
!  (front end wlSemiLeptonicOpacityModule2D) and the two-phase driver of
!  the momentum method (wlSemiLeptonicOpacityIntegrationModule2D) read
!  these tables; the integration kernels themselves are stateless
!  (Opacity_CC_2D_<Method>_Nodes take the rules as arguments).
!
!  Usage (host, before any offloaded region that calls Opacity_CC_2D):
!    CALL InitOpacity_CC_2D( nE [, Method_Option] )
!    ... Opacity_CC_2D( ... ) on host or device ...
!    CALL FinalizeOpacity_CC_2D()
!  Re-initialization (e.g. a different nE or method) is allowed between
!  kernels. The tables are read-only during the computation, so they are
!  safe for OpenMP threads and GPU threads alike; they are fixed-size
!  (nQuadMax nodes), which avoids allocatable declare-target data.
!
!  Integration methods (IntegrationMethods/):
!    iEmAb2D_Energy   = 1  ENERGY, the "standard" GSI-like method
!                          (outer and inner rule: gauleg(nE))
!    iEmAb2D_Momentum = 2  MOMENTUM (default): outer gauleg(nE), inner
!                          Inner_Rule(nE) = MAX(8,MIN(12,nE/2)) nodes
!    iEmAb2D_Tail     = 3  TAIL: gauleg(nE) for all pieces
!  The default of Method_Option, iEmAb2D_Default, is MOMENTUM unless the
!  library is built with EMAB2D_METHOD=ENERGY|TAIL (-DWEAKLIB_EMAB2D_*).
!=======================================================================
MODULE wlSemiLeptonicNodes2D

  USE wlKindModule, ONLY: dp
  USE wlSemiLeptonicPhysicsModule2D, ONLY: gauleg, nQuadMax
  USE wlSemiLeptonicIntegration2D_Momentum, ONLY: Inner_Rule, nE3Seg

  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER, PUBLIC :: iEmAb2D_Energy   = 1   ! "Standard" (GSI-like)
  INTEGER, PARAMETER, PUBLIC :: iEmAb2D_Momentum = 2
  INTEGER, PARAMETER, PUBLIC :: iEmAb2D_Tail     = 3
#if defined(WEAKLIB_EMAB2D_ENERGY)
  INTEGER, PARAMETER, PUBLIC :: iEmAb2D_Default  = iEmAb2D_Energy
#elif defined(WEAKLIB_EMAB2D_TAIL)
  INTEGER, PARAMETER, PUBLIC :: iEmAb2D_Default  = iEmAb2D_Tail
#else
  INTEGER, PARAMETER, PUBLIC :: iEmAb2D_Default  = iEmAb2D_Momentum
#endif

  ! ---- quadrature state (set by InitOpacity_CC_2D, read-only afterwards)
  LOGICAL , PUBLIC :: Initialized_CC2D = .FALSE.
  INTEGER , PUBLIC :: iMethod_CC2D     = iEmAb2D_Default
  INTEGER , PUBLIC :: nE_CC2D          = 0          ! resolution passed to Init
  ! outer rule (all methods)
  INTEGER , PUBLIC :: nE2_CC2D = 0
  REAL(DP), PUBLIC :: x2_CC2D(nQuadMax) = 0.0d0, w2_CC2D(nQuadMax) = 0.0d0
  ! inner rule of the selected method (ENERGY: gauleg(nE); MOMENTUM: Inner_Rule(nE))
  INTEGER , PUBLIC :: nE3_CC2D = 0
  REAL(DP), PUBLIC :: x3_CC2D(nQuadMax) = 0.0d0, w3_CC2D(nQuadMax) = 0.0d0
  ! inner rule of the MOMENTUM method, always set (two-phase driver)
  INTEGER , PUBLIC :: nE3M_CC2D = 0
  REAL(DP), PUBLIC :: x3M_CC2D(nE3Seg) = 0.0d0, w3M_CC2D(nE3Seg) = 0.0d0

#if defined(WEAKLIB_OMP_OL)
  !$OMP DECLARE TARGET( Initialized_CC2D, iMethod_CC2D, nE_CC2D,  &
  !$OMP                 nE2_CC2D, x2_CC2D, w2_CC2D,               &
  !$OMP                 nE3_CC2D, x3_CC2D, w3_CC2D,               &
  !$OMP                 nE3M_CC2D, x3M_CC2D, w3M_CC2D )
#elif defined(WEAKLIB_OACC)
  !$ACC DECLARE CREATE( Initialized_CC2D, iMethod_CC2D, nE_CC2D,  &
  !$ACC                 nE2_CC2D, x2_CC2D, w2_CC2D,               &
  !$ACC                 nE3_CC2D, x3_CC2D, w3_CC2D,               &
  !$ACC                 nE3M_CC2D, x3M_CC2D, w3M_CC2D )
#endif

  PUBLIC :: InitOpacity_CC_2D
  PUBLIC :: FinalizeOpacity_CC_2D

CONTAINS

!-----------------------------------------------------------------------
!  InitOpacity_CC_2D: compute the quadrature rules for resolution nE
!  (clamped to [1, nQuadMax]) and integration method Method_Option
!  (default iEmAb2D_Default; invalid values fall back to the default),
!  then copy them to the device. Host only.
!-----------------------------------------------------------------------
SUBROUTINE InitOpacity_CC_2D( nE, Method_Option )

  INTEGER, INTENT(IN)           :: nE
  INTEGER, INTENT(IN), OPTIONAL :: Method_Option

  INTEGER  :: nEuse
  REAL(DP) :: x3(nE3Seg), w3(nE3Seg)

  iMethod_CC2D = iEmAb2D_Default
  IF( PRESENT( Method_Option ) ) THEN
    IF( Method_Option == iEmAb2D_Energy .OR. Method_Option == iEmAb2D_Momentum &
        .OR. Method_Option == iEmAb2D_Tail ) iMethod_CC2D = Method_Option
  END IF

  nEuse   = MAX( 1, MIN( nE, nQuadMax ) )
  nE_CC2D = nEuse

  x2_CC2D = 0.0d0
  w2_CC2D = 0.0d0
  nE2_CC2D = nEuse
  CALL gauleg( 0.d0, 1.d0, x2_CC2D, w2_CC2D, nEuse )

  ! Momentum inner rule (also used by the two-phase driver)
  CALL Inner_Rule( nEuse, nE3M_CC2D, x3, w3 )
  x3M_CC2D = x3
  w3M_CC2D = w3

  x3_CC2D = 0.0d0
  w3_CC2D = 0.0d0
  SELECT CASE( iMethod_CC2D )
  CASE( iEmAb2D_Energy )
    nE3_CC2D = nEuse
    x3_CC2D  = x2_CC2D
    w3_CC2D  = w2_CC2D
  CASE( iEmAb2D_Momentum )
    nE3_CC2D = nE3M_CC2D
    x3_CC2D(1:nE3Seg) = x3M_CC2D
    w3_CC2D(1:nE3Seg) = w3M_CC2D
  CASE DEFAULT   ! Tail: no separate inner rule
    nE3_CC2D = 0
  END SELECT

  Initialized_CC2D = .TRUE.

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE TO( Initialized_CC2D, iMethod_CC2D, nE_CC2D,  &
  !$OMP                   nE2_CC2D, x2_CC2D, w2_CC2D,               &
  !$OMP                   nE3_CC2D, x3_CC2D, w3_CC2D,               &
  !$OMP                   nE3M_CC2D, x3M_CC2D, w3M_CC2D )
#elif defined(WEAKLIB_OACC)
  !$ACC UPDATE DEVICE( Initialized_CC2D, iMethod_CC2D, nE_CC2D,  &
  !$ACC                nE2_CC2D, x2_CC2D, w2_CC2D,               &
  !$ACC                nE3_CC2D, x3_CC2D, w3_CC2D,               &
  !$ACC                nE3M_CC2D, x3M_CC2D, w3M_CC2D )
#endif

END SUBROUTINE InitOpacity_CC_2D

!-----------------------------------------------------------------------
!  FinalizeOpacity_CC_2D: mark the tables as unset (Opacity_CC_2D then
!  returns 0). The arrays are static, so nothing is freed. Host only.
!-----------------------------------------------------------------------
SUBROUTINE FinalizeOpacity_CC_2D()

  Initialized_CC2D = .FALSE.
  nE_CC2D  = 0
  nE2_CC2D = 0
  nE3_CC2D = 0
  nE3M_CC2D = 0
  x3M_CC2D  = 0.0d0
  w3M_CC2D  = 0.0d0

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE TO( Initialized_CC2D, nE_CC2D, nE2_CC2D, nE3_CC2D, &
  !$OMP                   nE3M_CC2D, x3M_CC2D, w3M_CC2D )
#elif defined(WEAKLIB_OACC)
  !$ACC UPDATE DEVICE( Initialized_CC2D, nE_CC2D, nE2_CC2D, nE3_CC2D, &
  !$ACC                nE3M_CC2D, x3M_CC2D, w3M_CC2D )
#endif

END SUBROUTINE FinalizeOpacity_CC_2D

END MODULE wlSemiLeptonicNodes2D
