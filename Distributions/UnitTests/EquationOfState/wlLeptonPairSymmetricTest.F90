PROGRAM wlLeptonPairSymmetricTest

  USE, INTRINSIC :: IEEE_ARITHMETIC, ONLY: IEEE_IS_FINITE
  USE wlKindModule, ONLY: DP
  USE wlIOModuleHDF, ONLY: InitializeHDF, FinalizeHDF
  USE wlHelmIOModuleHDF, ONLY: ReadHelmholtzTableHDF
  USE wlLeptonEOSTableModule, ONLY: HelmTableType
  USE wlLeptonPhotonGasEOS, ONLY: &
    LeptonGasType, LeptonGasEOS, PairSymmetricEtaMin, PairSymmetricEtaMax, &
    ComputePairDegeneracy

  IMPLICIT NONE

  INTEGER, PARAMETER :: nRho = 4, nTemp = 5, nY = 97
  REAL(DP), PARAMETER :: RhoGrid(nRho) &
    = [ 1.0e7_DP, 1.0e10_DP, 1.0e12_DP, 1.0e14_DP ]
  REAL(DP), PARAMETER :: TempGrid(nTemp) &
    = [ 1.0e9_DP, 1.0e10_DP, 1.0e11_DP, 3.0e11_DP, 1.0e12_DP ]
  REAL(DP), PARAMETER :: SymmetryTolerance = 2.0e-12_DP
  REAL(DP), PARAMETER :: TransitionTolerance = 5.0e-5_DP
  REAL(DP), PARAMETER :: DerivativeTolerance = 2.0e-4_DP

  CHARACTER(256) :: TableName
  INTEGER :: iRho, iT, iY, nDerivative, nFailed, nFinite, nMonotone, nPositive
  INTEGER :: nSign, nSymmetry, nTransition, nTransitions, nZero
  REAL(DP) :: MaxTransitionJump, MuPrevious, RelativeJump, RelativeJumpInner
  REAL(DP) :: Rho, Temperature, Y
  TYPE(HelmTableType) :: MuonTable
  TYPE(LeptonGasType) :: MinusState, PlusState, ZeroState

  IF( COMMAND_ARGUMENT_COUNT() /= 1 )THEN
    ERROR STOP 'usage: wlLeptonPairSymmetricTest EOS_TABLE'
  END IF
  CALL GET_COMMAND_ARGUMENT( 1, TableName )

  CALL InitializeHDF
  CALL ReadHelmholtzTableHDF &
         ( MuonTable, TRIM(TableName), 'HelmTableMuons' )
  CALL FinalizeHDF

  nFailed      = 0
  nDerivative  = 0
  nFinite      = 0
  nMonotone    = 0
  nPositive    = 0
  nSign        = 0
  nSymmetry    = 0
  nTransition  = 0
  nTransitions = 0
  nZero        = 0
  MaxTransitionJump = 0.0_DP

  DO iT = 1, nTemp
    Temperature = TempGrid(iT)
    DO iRho = 1, nRho
      Rho = RhoGrid(iRho)

      CALL Evaluate( MuonTable, Rho, Temperature, 0.0_DP, ZeroState )
      IF( ZeroState % mu /= 0.0_DP ) nZero = nZero + 1
      IF( .NOT. StateIsFinite(ZeroState) ) nFinite = nFinite + 1
      IF( ZeroState % p < 0.0_DP .OR. ZeroState % e < 0.0_DP &
          .OR. ZeroState % s < 0.0_DP ) nPositive = nPositive + 1

      MuPrevious = 0.0_DP
      DO iY = 1, nY
        Y = 10.0_DP**( -30.0_DP &
                       + REAL(iY-1,DP) * 29.0_DP / REAL(nY-1,DP) )
        CALL Evaluate( MuonTable, Rho, Temperature,  Y, PlusState  )
        CALL Evaluate( MuonTable, Rho, Temperature, -Y, MinusState )

        IF( .NOT. StateIsFinite(PlusState) .OR. &
            .NOT. StateIsFinite(MinusState) ) nFinite = nFinite + 1
        IF( PlusState % mu <= 0.0_DP .OR. &
            MinusState % mu >= 0.0_DP ) nSign = nSign + 1
        IF( PlusState % mu < MuPrevious ) nMonotone = nMonotone + 1
        MuPrevious = PlusState % mu

        CALL CheckOddEvenSymmetry &
               ( PlusState, MinusState, SymmetryTolerance, nSymmetry )
      END DO

      CALL CheckTransition &
             ( MuonTable, Rho, Temperature, RelativeJump, nFinite )
      CALL CheckEtaBoundary &
             ( MuonTable, Rho, Temperature, PairSymmetricEtaMin, &
               RelativeJumpInner, nFinite )
      RelativeJump = MAX( RelativeJump, RelativeJumpInner )
      IF( RelativeJump >= 0.0_DP )THEN
        nTransitions = nTransitions + 1
        MaxTransitionJump = MAX( MaxTransitionJump, RelativeJump )
        IF( RelativeJump > TransitionTolerance )THEN
          nTransition = nTransition + 1
          WRITE(*,'(A,3ES12.4)') 'transition discontinuity: ', &
            Rho, Temperature, RelativeJump
        END IF
      END IF

      CALL CheckDerivatives &
             ( MuonTable, Rho, Temperature, DerivativeTolerance, &
               nDerivative )
    END DO
  END DO

  nFailed = nDerivative + nFinite + nMonotone + nPositive + nSign + nSymmetry &
            + nTransition + nZero
  WRITE(*,'(A,2ES12.4)') 'pair-symmetric eta transition: ', &
                          PairSymmetricEtaMin, PairSymmetricEtaMax
  WRITE(*,'(A,I0)') 'resolved transitions tested: ', nTransitions
  WRITE(*,'(A,ES12.4)') 'maximum transition relative jump: ', &
                         MaxTransitionJump
  WRITE(*,'(A,8(I0,1X))') &
    'derivative finite monotone positive sign symmetry transition zero failures: ', &
    nDerivative, nFinite, nMonotone, nPositive, nSign, nSymmetry, &
    nTransition, nZero
  WRITE(*,'(A,I0)') 'failed checks: ', nFailed
  IF( nFailed /= 0 ) ERROR STOP 'pair-symmetric EOS test failed'

CONTAINS

  SUBROUTINE Evaluate( Table, Rho, Temperature, Y, State )

    TYPE(HelmTableType), INTENT(IN) :: Table
    REAL(DP), INTENT(IN) :: Rho, Temperature, Y
    TYPE(LeptonGasType), INTENT(OUT) :: State

    State % rho = Rho
    State % T   = Temperature
    State % yL  = Y
    CALL LeptonGasEOS( Table, State )

  END SUBROUTINE Evaluate

  LOGICAL FUNCTION StateIsFinite( State )

    TYPE(LeptonGasType), INTENT(IN) :: State

    StateIsFinite &
      =     IEEE_IS_FINITE(State % p)    &
        .AND. IEEE_IS_FINITE(State % dpdT) &
        .AND. IEEE_IS_FINITE(State % dpdr) &
        .AND. IEEE_IS_FINITE(State % e)    &
        .AND. IEEE_IS_FINITE(State % dedT) &
        .AND. IEEE_IS_FINITE(State % dedr) &
        .AND. IEEE_IS_FINITE(State % s)    &
        .AND. IEEE_IS_FINITE(State % dsdT) &
        .AND. IEEE_IS_FINITE(State % dsdr) &
        .AND. IEEE_IS_FINITE(State % mu)

  END FUNCTION StateIsFinite

  SUBROUTINE CheckOddEvenSymmetry &
    ( PlusState, MinusState, Tolerance, nFailed )

    TYPE(LeptonGasType), INTENT(IN) :: PlusState, MinusState
    REAL(DP), INTENT(IN) :: Tolerance
    INTEGER, INTENT(INOUT) :: nFailed

    CALL CheckPair( PlusState % mu,  -MinusState % mu, Tolerance, nFailed )
    CALL CheckPair( PlusState % p,     MinusState % p, Tolerance, nFailed )
    CALL CheckPair( PlusState % dpdT,  MinusState % dpdT, Tolerance, nFailed )
    CALL CheckPair( PlusState % dpdr,  MinusState % dpdr, Tolerance, nFailed )
    CALL CheckPair( PlusState % e,     MinusState % e, Tolerance, nFailed )
    CALL CheckPair( PlusState % dedT,  MinusState % dedT, Tolerance, nFailed )
    CALL CheckPair( PlusState % dedr,  MinusState % dedr, Tolerance, nFailed )
    CALL CheckPair( PlusState % s,     MinusState % s, Tolerance, nFailed )
    CALL CheckPair( PlusState % dsdT,  MinusState % dsdT, Tolerance, nFailed )
    CALL CheckPair( PlusState % dsdr,  MinusState % dsdr, Tolerance, nFailed )

  END SUBROUTINE CheckOddEvenSymmetry

  SUBROUTINE CheckPair( A, B, Tolerance, nFailed )

    REAL(DP), INTENT(IN) :: A, B, Tolerance
    INTEGER, INTENT(INOUT) :: nFailed
    REAL(DP) :: Scale

    Scale = MAX( ABS(A), ABS(B), TINY(1.0_DP) )
    IF( ABS(A-B) > Tolerance * Scale ) nFailed = nFailed + 1

  END SUBROUTINE CheckPair

  SUBROUTINE CheckTransition &
    ( Table, Rho, Temperature, RelativeJump, nFailed )

    TYPE(HelmTableType), INTENT(IN) :: Table
    REAL(DP), INTENT(IN) :: Rho, Temperature
    REAL(DP), INTENT(OUT) :: RelativeJump
    INTEGER, INTENT(INOUT) :: nFailed

    INTEGER :: i
    LOGICAL :: HaveTransition
    REAL(DP) :: Eta, Jump, YFloor, YHi, YLo, YMid
    TYPE(LeptonGasType) :: HiState, LoState, MidState

    YLo = 0.0_DP
    YHi = 0.5_DP
    CALL ComputePairDegeneracy( Table, Temperature, Rho, YHi, Eta )
    IF( Eta <= PairSymmetricEtaMax )THEN
      RelativeJump = -1.0_DP
      RETURN
    END IF

    DO i = 1, 80
      YMid = 0.5_DP * ( YLo + YHi )
      CALL ComputePairDegeneracy &
             ( Table, Temperature, Rho, YMid, Eta )
      IF( Eta <= PairSymmetricEtaMax )THEN
        YLo = YMid
      ELSE
        YHi = YMid
      END IF
    END DO

    RelativeJump = -1.0_DP
    HaveTransition = .FALSE.
    IF( Rho*YMid >= Table % mindens )THEN
      CALL Evaluate( Table, Rho, Temperature, &
                     YLo*(1.0_DP-1.0e-7_DP), LoState )
      CALL Evaluate( Table, Rho, Temperature, &
                     YHi*(1.0_DP+1.0e-7_DP), HiState )
      RelativeJump = StateRelativeDifference( HiState, LoState )
      HaveTransition = .TRUE.
    END IF

    YFloor = Table % mindens / Rho
    IF( YFloor > 0.0_DP .AND. YFloor < 0.5_DP )THEN
      CALL ComputePairDegeneracy &
             ( Table, Temperature, Rho, YFloor, Eta )
      IF( Eta > PairSymmetricEtaMax )THEN
        CALL Evaluate( Table, Rho, Temperature, &
                       YFloor*(1.0_DP-1.0e-7_DP), LoState )
        CALL Evaluate( Table, Rho, Temperature, &
                       YFloor*(1.0_DP+1.0e-7_DP), HiState )
        Jump = StateRelativeDifference( HiState, LoState )
        RelativeJump = MAX( RelativeJump, Jump )
        HaveTransition = .TRUE.
      END IF
    END IF
    IF( HaveTransition )THEN
      IF( .NOT. StateIsFinite(LoState) .OR. &
          .NOT. StateIsFinite(HiState) ) nFailed = nFailed + 1
    END IF

  END SUBROUTINE CheckTransition

  SUBROUTINE CheckEtaBoundary &
    ( Table, Rho, Temperature, TargetEta, RelativeJump, nFailed )

    TYPE(HelmTableType), INTENT(IN) :: Table
    REAL(DP), INTENT(IN) :: Rho, Temperature, TargetEta
    REAL(DP), INTENT(OUT) :: RelativeJump
    INTEGER, INTENT(INOUT) :: nFailed

    INTEGER :: i
    REAL(DP) :: Eta, YHi, YLo, YMid
    TYPE(LeptonGasType) :: HiState, LoState

    YLo = 0.0_DP
    YHi = 0.5_DP
    CALL ComputePairDegeneracy( Table, Temperature, Rho, YHi, Eta )
    IF( Eta <= TargetEta )THEN
      RelativeJump = -1.0_DP
      RETURN
    END IF
    DO i = 1, 80
      YMid = 0.5_DP * ( YLo + YHi )
      CALL ComputePairDegeneracy( Table, Temperature, Rho, YMid, Eta )
      IF( Eta <= TargetEta )THEN
        YLo = YMid
      ELSE
        YHi = YMid
      END IF
    END DO

    RelativeJump = -1.0_DP
    IF( Rho*YMid >= Table % mindens )THEN
      CALL Evaluate( Table, Rho, Temperature, &
                     YLo*(1.0_DP-1.0e-7_DP), LoState )
      CALL Evaluate( Table, Rho, Temperature, &
                     YHi*(1.0_DP+1.0e-7_DP), HiState )
      RelativeJump = StateRelativeDifference( HiState, LoState )
      IF( .NOT. StateIsFinite(LoState) .OR. &
          .NOT. StateIsFinite(HiState) ) nFailed = nFailed + 1
    END IF

  END SUBROUTINE CheckEtaBoundary

  REAL(DP) FUNCTION StateRelativeDifference( A, B )

    TYPE(LeptonGasType), INTENT(IN) :: A, B

    StateRelativeDifference = 0.0_DP
    CALL IncludeDifference( A % p,    B % p,    StateRelativeDifference )
    CALL IncludeDifference( A % e,    B % e,    StateRelativeDifference )
    CALL IncludeDifference( A % s,    B % s,    StateRelativeDifference )
    CALL IncludeDifference( A % mu,   B % mu,   StateRelativeDifference )

  END FUNCTION StateRelativeDifference

  SUBROUTINE IncludeDifference( A, B, Maximum )

    REAL(DP), INTENT(IN) :: A, B
    REAL(DP), INTENT(INOUT) :: Maximum
    REAL(DP) :: Scale

    Scale = MAX( ABS(A), ABS(B), TINY(1.0_DP) )
    Maximum = MAX( Maximum, ABS(A-B)/Scale )

  END SUBROUTINE IncludeDifference

  SUBROUTINE CheckDerivatives &
    ( Table, Rho, Temperature, Tolerance, nFailed )

    TYPE(HelmTableType), INTENT(IN) :: Table
    REAL(DP), INTENT(IN) :: Rho, Temperature, Tolerance
    INTEGER, INTENT(INOUT) :: nFailed

    INTEGER :: i
    REAL(DP), PARAMETER :: h = 1.0e-4_DP
    REAL(DP) :: Eta, YHi, YLo, YMid
    TYPE(LeptonGasType) :: Center, Minus, Plus

    YLo = 0.0_DP
    YHi = 0.5_DP
    CALL ComputePairDegeneracy( Table, Temperature, Rho, YHi, Eta )
    IF( Eta <= 0.5_DP*PairSymmetricEtaMin ) RETURN

    DO i = 1, 80
      YMid = 0.5_DP * ( YLo + YHi )
      CALL ComputePairDegeneracy &
             ( Table, Temperature, Rho, YMid, Eta )
      IF( Eta <= 0.5_DP*PairSymmetricEtaMin )THEN
        YLo = YMid
      ELSE
        YHi = YMid
      END IF
    END DO
    YMid = 0.5_DP * ( YLo + YHi )

    CALL Evaluate( Table, Rho, Temperature, YMid, Center )
    CALL Evaluate( Table, Rho*(1.0_DP-h), Temperature, YMid, Minus )
    CALL Evaluate( Table, Rho*(1.0_DP+h), Temperature, YMid, Plus )
    CALL CheckDerivative &
           ( Rho*Center % dpdr, (Plus % p-Minus % p)/(2.0_DP*h), &
             Center % p, h, Tolerance, 'rho dp/drho', nFailed )
    CALL CheckDerivative &
           ( Rho*Center % dedr, (Plus % e-Minus % e)/(2.0_DP*h), &
             Center % e, h, Tolerance, 'rho de/drho', nFailed )
    CALL CheckDerivative &
           ( Rho*Center % dsdr, (Plus % s-Minus % s)/(2.0_DP*h), &
             Center % s, h, Tolerance, 'rho ds/drho', nFailed )

    CALL Evaluate( Table, Rho, Temperature*(1.0_DP-h), YMid, Minus )
    CALL Evaluate( Table, Rho, Temperature*(1.0_DP+h), YMid, Plus )
    CALL CheckDerivative &
           ( Temperature*Center % dpdT, &
             (Plus % p-Minus % p)/(2.0_DP*h), &
             Center % p, h, Tolerance, 'T dp/dT', nFailed )
    CALL CheckDerivative &
           ( Temperature*Center % dedT, &
             (Plus % e-Minus % e)/(2.0_DP*h), &
             Center % e, h, Tolerance, 'T de/dT', nFailed )
    CALL CheckDerivative &
           ( Temperature*Center % dsdT, &
             (Plus % s-Minus % s)/(2.0_DP*h), &
             Center % s, h, Tolerance, 'T ds/dT', nFailed )

  END SUBROUTINE CheckDerivatives

  SUBROUTINE CheckDerivative &
    ( Analytic, Difference, Reference, h, Tolerance, Label, nFailed )

    REAL(DP), INTENT(IN) :: Analytic, Difference, Reference, h, Tolerance
    CHARACTER(*), INTENT(IN) :: Label
    INTEGER, INTENT(INOUT) :: nFailed
    REAL(DP) :: RoundoffFloor, Scale

    Scale = MAX( ABS(Analytic), ABS(Difference), TINY(1.0_DP) )
    RoundoffFloor = 32.0_DP * EPSILON(1.0_DP) * ABS(Reference) / h
    IF( ABS(Analytic-Difference) > &
        Tolerance*Scale + RoundoffFloor )THEN
      nFailed = nFailed + 1
      WRITE(*,'(A,1X,A,3ES14.6)') 'derivative mismatch:', Label, &
        Analytic, Difference, ABS(Analytic-Difference)/Scale
    END IF

  END SUBROUTINE CheckDerivative

END PROGRAM wlLeptonPairSymmetricTest
