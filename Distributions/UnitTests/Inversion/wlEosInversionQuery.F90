PROGRAM wlEosInversionQuery

  USE wlKindModule, ONLY: &
    dp
  USE wlEquationOfStateTableModule, ONLY: &
    EquationOfStateTableType
  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF
  USE wlEOSIOModuleHDF, ONLY: &
    ReadEquationOfStateTableHDF
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable_3D_Custom, &
    LogInterpolateSingleVariable_3D_Custom_Point
  USE wlEOSInversionModule, ONLY: &
    InitializeEOSInversion, &
    ComputeTemperatureWith_DEY_Many, &
    ComputeTemperatureWith_DEY_Single_Guess, &
    ComputeTemperatureWith_DEY_Single_NoGuess, &
    ComputeTemperatureWith_DSY_Many, &
    ComputeTemperatureWith_DSY_Single_Guess, &
    ComputeTemperatureWith_DSY_Single_NoGuess, &
    ComputeTemperatureWith_DPY_Many, &
    ComputeTemperatureWith_DPY_Single_Guess, &
    ComputeTemperatureWith_DPY_Single_NoGuess, &
    DescribeEOSInversionError

  IMPLICIT NONE

  REAL(dp), PARAMETER :: D_0 = 1.032E12
  REAL(dp), PARAMETER :: E_0 = 2.3173722677362975E19
  REAL(dp), PARAMETER :: Y_0 = 0.1347

  INTEGER :: &
    n_rndm, &
    iP, &
    iMaxError
  INTEGER, PARAMETER :: &
    nPoints = 2**16, &
    iD = 1, iT = 2, iY = 3
  INTEGER, DIMENSION(nPoints) :: &
    Error_E, &
    Error_P, &
    Error_S
  REAL(dp) :: &
    tBegin, &
    tEnd, &
    tCPU, &
    tGPU, &
    T_Guess, &
    Amp, &
    E_T, S_T, P_T
  REAL(dp), DIMENSION(nPoints) :: &
    D, T, Y, E, P, S, T_E, T_P, T_S, dT, &
    rndm_D, &
    rndm_T, &
    rndm_Y, &
    Error
  TYPE(EquationOfStateTableType) :: &
    EOS
  REAL(DP), DIMENSION(:), ALLOCATABLE :: &
    Ds_T, Ts_T, Ys_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: &
    Ps_T, Es_T, Ss_T
  REAL(DP) :: &
    OS_P, OS_E, OS_S
  INTEGER :: &
    iD_T, iT_T, iY_T, iP_T, iS_T, iE_T

#if defined(WEAKLIB_OMP_OL)
#elif defined(WEAKLIB_OACC)
  !$ACC INIT
#endif

  CALL InitializeHDF( )
  CALL ReadEquationOfStateTableHDF( EOS, "EquationOfStateTable.h5" )
  CALL FinalizeHDF( )

  iD_T = EOS % TS % Indices % iRho
  iT_T = EOS % TS % Indices % iT
  iY_T = EOS % TS % Indices % iYe

  ALLOCATE( Ds_T(EOS % TS % nPoints(iD_T)) )
  Ds_T = EOS % TS % States(iD_T) % Values

  ALLOCATE( Ts_T(EOS % TS % nPoints(iT_T)) )
  Ts_T = EOS % TS % States(iT_T) % Values

  ALLOCATE( Ys_T(EOS % TS % nPoints(iY_T)) )
  Ys_T = EOS % TS % States(iY_T) % Values

  iP_T = EOS % DV % Indices % iPressure
  iS_T = EOS % DV % Indices % iEntropyPerBaryon
  iE_T = EOS % DV % Indices % iInternalEnergyDensity

  OS_P = EOS % DV % Offsets(iP_T)
  OS_S = EOS % DV % Offsets(iS_T)
  OS_E = EOS % DV % Offsets(iE_T)

  ALLOCATE &
    ( Ps_T(1:EOS % DV % nPoints(1), &
           1:EOS % DV % nPoints(2), &
           1:EOS % DV % nPoints(3)) )
  ALLOCATE &
    ( Ss_T(1:EOS % DV % nPoints(1), &
           1:EOS % DV % nPoints(2), &
           1:EOS % DV % nPoints(3)) )
  ALLOCATE &
    ( Es_T(1:EOS % DV % nPoints(1), &
           1:EOS % DV % nPoints(2), &
           1:EOS % DV % nPoints(3)) )

  Ps_T = EOS % DV % Variables(iP_T ) % Values
  Ss_T = EOS % DV % Variables(iS_T ) % Values
  Es_T = EOS % DV % Variables(iE_T ) % Values

  CALL InitializeEOSInversion &
         ( Ds_T, Ts_T, Ys_T, &
           10.0d0**( Es_T ) - OS_E, &
           10.0d0**( Ps_T ) - OS_P, &
           10.0d0**( Ss_T ) - OS_S, &
           Verbose_Option = .TRUE. )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET ENTER DATA &
  !$OMP MAP( to: Ds_T, Ts_T, Ys_T, Es_T, Ps_T, Ss_T, OS_E, OS_P, OS_S ) &
  !$OMP MAP( alloc: D, T, Y, P, E, S, T_P, T_E, T_S, Error_P, Error_E, Error_S )
#elif defined (WEAKLIB_OACC)
  !$ACC ENTER DATA &
  !$ACC COPYIN( Ds_T, Ts_T, Ys_T, Es_T, Ps_T, Ss_T, OS_E, OS_P, OS_S ) &
  !$ACC CREATE( D, T, Y, P, E, S, T_P, T_E, T_S, Error_P, Error_E, Error_S )
#endif

  WRITE(*,*)
  WRITE(*,'(A4,A10,I10.10)') '', 'nPoints = ', nPoints

  ! --- Random Sampling of EOS Table ---

  n_rndm = nPoints

  ! --- Initialize Density Points ---

  D = D_0

  PRINT*, "Min/Max D = ", MINVAL( D ), MAXVAL( D )

  ! --- Initialize Energy Points ---

  E = E_0

  PRINT*, "Min/Max E = ", MINVAL( E ), MAXVAL( E )

  ! --- Initialize Electron Fraction Points ---

  Y = Y_0

  PRINT*, "Min/Max Y = ", MINVAL( Y ), MAXVAL( Y )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE TO( D, Y, E )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE DEVICE( D, Y, E )
#endif

  ! -------------------------------------------------------------------

  T_E = 0.0_dp
  Error_E = 0

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE TO( T_E, Error_E )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE DEVICE( T_E, Error_E )
#endif

  CALL CPU_TIME( tBegin )

  CALL ComputeTemperatureWith_DEY_Many &
         ( D, E, Y, Ds_T, Ts_T, Ys_T, Es_T, OS_E, T_E, &
           UseInitialGuess_Option = .FALSE., &
           Error_Option = Error_E )

  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  CALL CPU_TIME( tBegin )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined (WEAKLIB_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR &
  !$ACC PRESENT( D, E, Y, Ds_T, Ts_T, Ys_T, Es_T, OS_E, T_E, Error_E )
#endif
  DO iP = 1, nPoints
    CALL ComputeTemperatureWith_DEY_Single_NoGuess &
           ( D(iP), E(iP), Y(iP), Ds_T, Ts_T, Ys_T, Es_T, OS_E, T_E(iP), &
             Error_E(iP) )
  END DO

  CALL CPU_TIME( tEnd )
  tGPU = tEnd - tBegin

  PRINT*
  PRINT*, "ComputeTemperatureWith_DEY (No Guess):"
  PRINT*, "CPU_TIME = ", tCPU, "GPU_TIME = ", tGPU

  DO iP = 1, nPoints
    IF( Error_E(iP) .NE. 0 ) CALL DescribeEOSInversionError( Error_E(iP) )
  END DO
  PRINT*, "CPU Min/Max T = ", MINVAL( T_E ), MAXVAL( T_E )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE FROM( T_E, Error_E )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE HOST( T_E, Error_E )
#endif
  DO iP = 1, nPoints
    IF( Error_E(iP) .NE. 0 ) CALL DescribeEOSInversionError( Error_E(iP) )
  END DO
  PRINT*, "GPU Min/Max T = ", MINVAL( T_E ), MAXVAL( T_E )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET EXIT DATA &
  !$OMP MAP( release: D, E, P, S, Y, Ds_T, Ts_T, Ys_T, Es_T, Ps_T, Ss_T, OS_E, OS_P, OS_S, T_E, Error_E )
#elif defined (WEAKLIB_OACC)
  !$ACC EXIT DATA &
  !$ACC DELETE( D, E, P, S, Y, Ds_T, Ts_T, Ys_T, Es_T, Ps_T, Ss_T, OS_E, OS_P, OS_S, T_E, Error_E )
#endif

  DEALLOCATE( Ds_T, Ts_T, Ys_T, Ps_T, Ss_T, Es_T )

#if defined(WEAKLIB_OMP_OL)
#elif defined(WEAKLIB_OACC)
  !$ACC SHUTDOWN
#endif

END PROGRAM wlEosInversionQuery
