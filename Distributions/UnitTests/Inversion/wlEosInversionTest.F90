PROGRAM wlEosInversionTest

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

  INCLUDE 'mpif.h'

  INTEGER :: &
    n_rndm, &
    iP, &
    iMaxError, &
    ierr
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

  CALL MPI_INIT( ierr )

#if defined(WEAKLIB_OMP_OL)
#elif defined(WEAKLIB_OACC)
  !!$ACC INIT
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

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_D )

  associate &
    ( minD => EOS % TS % MinValues(iD), &
      maxD => EOS % TS % MaxValues(iD) )

  D(:) = 10**( LOG10(minD) + ( LOG10(maxD) - LOG10(minD) ) * rndm_D )

  PRINT*, "Min/Max D = ", MINVAL( D ), MAXVAL( D )
  PRINT*, "            ", minD, maxD

  end associate

  ! --- Initialize Temperature Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_T )

  associate &
    ( minT => EOS % TS % MinValues(iT), &
      maxT => EOS % TS % MaxValues(iT) )

  T(:) = 10**( LOG10(minT) + ( LOG10(maxT) - LOG10(minT) ) * rndm_T )

  PRINT*, "Min/Max T = ", MINVAL( T ), MAXVAL( T )
  PRINT*, "            ", minT, maxT

  end associate

  ! --- Initialize Electron Fraction Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_Y )

  associate &
    ( minY => EOS % TS % MinValues(iY), &
      maxY => EOS % TS % MaxValues(iY) )

  Y(:) = minY + ( maxY - minY ) * rndm_Y

  PRINT*, "Min/Max Y = ", MINVAL( Y ), MAXVAL( Y )
  PRINT*, "            ", minY, maxY

  end associate

  ! --- Compute Internal Energy, Pressure, and Entropy ---

  DO iP = 1, nPoints

    CALL LogInterpolateSingleVariable_3D_Custom_Point &
           ( D(iP), T(iP), Y(iP), Ds_T, Ts_T, Ys_T, OS_E, Es_T, E(iP) )

    CALL LogInterpolateSingleVariable_3D_Custom_Point &
           ( D(iP), T(iP), Y(iP), Ds_T, Ts_T, Ys_T, OS_P, Ps_T, P(iP) )

    CALL LogInterpolateSingleVariable_3D_Custom_Point &
           ( D(iP), T(iP), Y(iP), Ds_T, Ts_T, Ys_T, OS_S, Ss_T, S(iP) )

  END DO

  WRITE(*,*)
  WRITE(*,*) "Min/Max E = ", MINVAL( E ), MAXVAL( E )
  WRITE(*,*) "Min/Max P = ", MINVAL( P ), MAXVAL( P )
  WRITE(*,*) "Min/Max S = ", MINVAL( S ), MAXVAL( S )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE TO( D, T, Y, P, E, S )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE DEVICE( D, T, Y, P, E, S )
#endif

  ! -------------------------------------------------------------------
  ! --- Recover Temperature from Internal Energy ----------------------
  ! -------------------------------------------------------------------

  Amp = 0.001_dp ! --- Perturbation Amplitude

  PRINT*
  PRINT*, "Perturbation Amplitude: ", Amp

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_T )

  dT = Amp * 2.0_dp * ( rndm_T - 0.5_dp ) * T

  associate &
    ( minT => 1.0001 * EOS % TS % MinValues(iT), &
      maxT => 0.9999 * EOS % TS % MaxValues(iT) )

  T_E = MIN( MAX( T + dT, minT ), maxT )

  end associate

  Error_E = 0

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE TO( T_E, Error_E )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE DEVICE( T_E, Error_E )
#endif

  CALL CPU_TIME( tBegin )

  CALL ComputeTemperatureWith_DEY_Many &
         ( D, E, Y, Ds_T, Ts_T, Ys_T, Es_T, OS_E, T_E, &
           UseInitialGuess_Option = .TRUE., &
           Error_Option = Error_E )

  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  CALL CPU_TIME( tBegin )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  !$OMP PRIVATE( T_Guess )
#elif defined (WEAKLIB_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR &
  !$ACC PRIVATE( T_Guess ) &
  !$ACC PRESENT( D, E, Y, Ds_T, Ts_T, Ys_T, Es_T, OS_E, T_E, Error_E )
#elif defined (WEAKLIB_OMP)
  !$OMP PARALLEL DO &
  !$OMP PRIVATE( T_Guess )
#endif
  DO iP = 1, nPoints
    T_Guess = T_E(iP)
    CALL ComputeTemperatureWith_DEY_Single_Guess &
           ( D(iP), E(iP), Y(iP), Ds_T, Ts_T, Ys_T, Es_T, OS_E, T_E(iP), T_Guess, &
             Error_Option = Error_E(iP) )
  END DO

  CALL CPU_TIME( tEnd )
  tGPU = tEnd - tBegin

  PRINT*
  PRINT*, "ComputeTemperatureWith_DEY (Good Guess):"
  PRINT*, "CPU_TIME = ", tCPU, "GPU_TIME = ", tGPU

  DO iP = 1, nPoints
    IF( Error_E(iP) .NE. 0 ) CALL DescribeEOSInversionError( Error_E(iP) )
  END DO
  Error = ABS( T - T_E ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_E(iMaxError), Y(iMaxError), Ds_T, Ts_T, Ys_T, OS_E, Es_T, E_T )
  PRINT*, "CPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "CPU  Error(E) = ", ABS( E(iMaxError) - E_T ) / E(iMaxError)
  PRINT*, "CPU MINVAL(T) = ", MINVAL( T_E )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE FROM( T_E, Error_E )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE HOST( T_E, Error_E )
#endif
  DO iP = 1, nPoints
    IF( Error_E(iP) .NE. 0 ) CALL DescribeEOSInversionError( Error_E(iP) )
  END DO
  Error = ABS( T - T_E ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_E(iMaxError), Y(iMaxError), Ds_T, Ts_T, Ys_T, OS_E, Es_T, E_T )
  PRINT*, "GPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "GPU  Error(E) = ", ABS( E(iMaxError) - E_T ) / E(iMaxError)
  PRINT*, "GPU MINVAL(T) = ", MINVAL( T_E )

  ! -------------------------------------------------------------------

  associate &
    ( minT => 1.0001 * EOS % TS % MinValues(iT), &
      maxT => 0.9999 * EOS % TS % MaxValues(iT) )

  Amp = 0.1_dp * ( maxT - minT ) ! --- Perturbation Amplitude

  PRINT*
  PRINT*, "Perturbation Amplitude: ", Amp

  dT = Amp * 2.0_dp * ( rndm_T - 0.5_dp ) * T

  T_E = MIN( MAX( T + dT, minT ), maxT )

  end associate

  Error_E = 0

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE TO( T_E, Error_E )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE DEVICE( T_E, Error_E )
#endif

  CALL CPU_TIME( tBegin )

  CALL ComputeTemperatureWith_DEY_Many &
         ( D, E, Y, Ds_T, Ts_T, Ys_T, Es_T, OS_E, T_E, &
           UseInitialGuess_Option = .TRUE., &
           Error_Option = Error_E )

  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  CALL CPU_TIME( tBegin )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  !$OMP PRIVATE( T_Guess )
#elif defined (WEAKLIB_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR &
  !$ACC PRIVATE( T_Guess ) &
  !$ACC PRESENT( D, E, Y, Ds_T, Ts_T, Ys_T, Es_T, OS_E, T_E, Error_E )
#elif defined (WEAKLIB_OMP)
  !$OMP PARALLEL DO &
  !$OMP PRIVATE( T_Guess )
#endif
  DO iP = 1, nPoints
    T_Guess = T_E(iP)
    CALL ComputeTemperatureWith_DEY_Single_Guess &
           ( D(iP), E(iP), Y(iP), Ds_T, Ts_T, Ys_T, Es_T, OS_E, T_E(iP), T_Guess, &
             Error_Option = Error_E(iP) )
  END DO

  CALL CPU_TIME( tEnd )
  tGPU = tEnd - tBegin

  PRINT*
  PRINT*, "ComputeTemperatureWith_DEY (Bad Guess):"
  PRINT*, "CPU_TIME = ", tCPU, "GPU_TIME = ", tGPU

  DO iP = 1, nPoints
    IF( Error_E(iP) .NE. 0 ) CALL DescribeEOSInversionError( Error_E(iP) )
  END DO
  Error = ABS( T - T_E ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_E(iMaxError), Y(iMaxError), Ds_T, Ts_T, Ys_T, OS_E, Es_T, E_T )
  PRINT*, "CPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "CPU  Error(E) = ", ABS( E(iMaxError) - E_T ) / E(iMaxError)
  PRINT*, "CPU MINVAL(T) = ", MINVAL( T_E )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE FROM( T_E, Error_E )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE HOST( T_E, Error_E )
#endif
  DO iP = 1, nPoints
    IF( Error_E(iP) .NE. 0 ) CALL DescribeEOSInversionError( Error_E(iP) )
  END DO
  Error = ABS( T - T_E ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_E(iMaxError), Y(iMaxError), Ds_T, Ts_T, Ys_T, OS_E, Es_T, E_T )
  PRINT*, "GPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "GPU  Error(E) = ", ABS( E(iMaxError) - E_T ) / E(iMaxError)
  PRINT*, "GPU MINVAL(T) = ", MINVAL( T_E )

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
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  !$OMP PRIVATE( T_Guess )
#elif defined (WEAKLIB_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR &
  !$ACC PRIVATE( T_Guess ) &
  !$ACC PRESENT( D, E, Y, Ds_T, Ts_T, Ys_T, Es_T, OS_E, T_E, Error_E )
#elif defined (WEAKLIB_OMP)
  !$OMP PARALLEL DO &
  !$OMP PRIVATE( T_Guess )
#endif
  DO iP = 1, nPoints
    T_Guess = T_E(iP)
    CALL ComputeTemperatureWith_DEY_Single_NoGuess &
           ( D(iP), E(iP), Y(iP), Ds_T, Ts_T, Ys_T, Es_T, OS_E, T_E(iP), &
             Error_Option = Error_E(iP) )
  END DO

  CALL CPU_TIME( tEnd )
  tGPU = tEnd - tBegin

  PRINT*
  PRINT*, "ComputeTemperatureWith_DEY (No Guess):"
  PRINT*, "CPU_TIME = ", tCPU, "GPU_TIME = ", tGPU

  DO iP = 1, nPoints
    IF( Error_E(iP) .NE. 0 ) CALL DescribeEOSInversionError( Error_E(iP) )
  END DO
  Error = ABS( T - T_E ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_E(iMaxError), Y(iMaxError), Ds_T, Ts_T, Ys_T, OS_E, Es_T, E_T )
  PRINT*, "CPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "CPU  Error(E) = ", ABS( E(iMaxError) - E_T ) / E(iMaxError)
  PRINT*, "CPU MINVAL(T) = ", MINVAL( T_E )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE FROM( T_E, Error_E )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE HOST( T_E, Error_E )
#endif
  DO iP = 1, nPoints
    IF( Error_E(iP) .NE. 0 ) CALL DescribeEOSInversionError( Error_E(iP) )
  END DO
  Error = ABS( T - T_E ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_E(iMaxError), Y(iMaxError), Ds_T, Ts_T, Ys_T, OS_E, Es_T, E_T )
  PRINT*, "GPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "GPU  Error(E) = ", ABS( E(iMaxError) - E_T ) / E(iMaxError)
  PRINT*, "GPU MINVAL(T) = ", MINVAL( T_E )

  ! -------------------------------------------------------------------

  !T_E = 0.0_dp

  !CALL CPU_TIME( tBegin )

  !CALL ComputeTempFromIntEnergy_Lookup &
  !       ( D, E, Y, Ds_T, Ts_T, Ys_T, [1,1,0], &
  !         Es_T, OS_E, T_E )

  !CALL CPU_TIME( tEnd )

  !Error = ABS( T - T_E ) / T
  !iMaxError = MAXLOC( Error, DIM=1 )

  !CALL LogInterpolateSingleVariable_3D_Custom_Point &
  !       ( D(iMaxError), T_E(iMaxError), Y(iMaxError), Ds_T, Ts_T, Ys_T, OS_E, Es_T, E_T )

  !PRINT*
  !PRINT*, "ComputeTempFromIntEnergy_Lookup:"
  !PRINT*, "CPU_TIME = ", tEnd - tBegin
  !PRINT*, MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  !PRINT*, "MINVAL(T) = ", MINVAL( T_E )

  !STOP

  ! -------------------------------------------------------------------
  ! --- Recover Temperature from Entropy Per Baryon -------------------
  ! -------------------------------------------------------------------

  Amp = 0.001_dp
  dT = Amp * 2.0_dp * ( rndm_T - 0.5_dp ) * T

  associate &
    ( minT => 1.0001 * EOS % TS % MinValues(iT), &
      maxT => 0.9999 * EOS % TS % MaxValues(iT) )

  T_S = MIN( MAX( T + dT, minT ), maxT )

  end associate

  Error_S = 0

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE TO( T_S, Error_S )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE DEVICE( T_S, Error_S )
#endif

  CALL CPU_TIME( tBegin )

  CALL ComputeTemperatureWith_DSY_Many &
         ( D, S, Y, Ds_T, Ts_T, Ys_T, Ss_T, OS_S, T_S, &
           UseInitialGuess_Option = .TRUE., &
           Error_Option = Error_S )

  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  CALL CPU_TIME( tBegin )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  !$OMP PRIVATE( T_Guess )
#elif defined (WEAKLIB_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR &
  !$ACC PRIVATE( T_Guess ) &
  !$ACC PRESENT( D, S, Y, Ds_T, Ts_T, Ys_T, Ss_T, OS_S, T_S, Error_S )
#elif defined (WEAKLIB_OMP)
  !$OMP PARALLEL DO &
  !$OMP PRIVATE( T_Guess )
#endif
  DO iP = 1, nPoints
    T_Guess = T_S(iP)
    CALL ComputeTemperatureWith_DSY_Single_Guess &
           ( D(iP), S(iP), Y(iP), Ds_T, Ts_T, Ys_T, Ss_T, OS_S, T_S(iP), T_Guess, &
             Error_Option = Error_S(iP) )
  END DO

  CALL CPU_TIME( tEnd )
  tGPU = tEnd - tBegin

  PRINT*
  PRINT*, "ComputeTemperatureWith_DSY (Good Guess):"
  PRINT*, "CPU_TIME = ", tCPU, "GPU_TIME = ", tGPU

  DO iP = 1, nPoints
    IF( Error_S(iP) .NE. 0 ) CALL DescribeEOSInversionError( Error_S(iP) )
  END DO
  Error = ABS( T - T_S ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_S(iMaxError), Y(iMaxError), Ds_T, Ts_T, Ys_T, OS_S, Ss_T, S_T )
  PRINT*, "CPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "CPU  Error(S) = ", ABS( S(iMaxError) - S_T ) / S(iMaxError)
  PRINT*, "CPU MINVAL(T) = ", MINVAL( T_S )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE FROM( T_S, Error_S )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE HOST( T_S, Error_S )
#endif
  DO iP = 1, nPoints
    IF( Error_S(iP) .NE. 0 ) CALL DescribeEOSInversionError( Error_S(iP) )
  END DO
  Error = ABS( T - T_S ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_S(iMaxError), Y(iMaxError), Ds_T, Ts_T, Ys_T, OS_S, Ss_T, S_T )
  PRINT*, "GPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "GPU  Error(S) = ", ABS( S(iMaxError) - S_T ) / S(iMaxError)
  PRINT*, "GPU MINVAL(T) = ", MINVAL( T_S )

  ! -------------------------------------------------------------------

  associate &
    ( minT => 1.0001 * EOS % TS % MinValues(iT), &
      maxT => 0.9999 * EOS % TS % MaxValues(iT) )

  Amp = 0.1_dp * ( maxT - minT )
  dT = Amp * 2.0_dp * ( rndm_T - 0.5_dp ) * T

  T_S = MIN( MAX( T + dT, minT ), maxT )

  end associate

  Error_S = 0

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE TO( T_S, Error_S )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE DEVICE( T_S, Error_S )
#endif

  CALL CPU_TIME( tBegin )

  CALL ComputeTemperatureWith_DSY_Many &
         ( D, S, Y, Ds_T, Ts_T, Ys_T, Ss_T, OS_S, T_S, &
           UseInitialGuess_Option = .TRUE., &
           Error_Option = Error_S )

  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  CALL CPU_TIME( tBegin )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  !$OMP PRIVATE( T_Guess )
#elif defined (WEAKLIB_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR &
  !$ACC PRIVATE( T_Guess ) &
  !$ACC PRESENT( D, S, Y, Ds_T, Ts_T, Ys_T, Ss_T, OS_S, T_S, Error_S )
#elif defined (WEAKLIB_OMP)
  !$OMP PARALLEL DO &
  !$OMP PRIVATE( T_Guess )
#endif
  DO iP = 1, nPoints
    T_Guess = T_S(iP)
    CALL ComputeTemperatureWith_DSY_Single_Guess &
           ( D(iP), S(iP), Y(iP), Ds_T, Ts_T, Ys_T, Ss_T, OS_S, T_S(iP), T_Guess, &
             Error_Option = Error_S(iP) )
  END DO

  CALL CPU_TIME( tEnd )
  tGPU = tEnd - tBegin

  PRINT*
  PRINT*, "ComputeTemperatureWith_DSY (Bad Guess):"
  PRINT*, "CPU_TIME = ", tCPU, "GPU_TIME = ", tGPU

  DO iP = 1, nPoints
    IF( Error_S(iP) .NE. 0 ) CALL DescribeEOSInversionError( Error_S(iP) )
  END DO
  Error = ABS( T - T_S ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_S(iMaxError), Y(iMaxError), Ds_T, Ts_T, Ys_T, OS_S, Ss_T, S_T )
  PRINT*, "CPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "CPU  Error(S) = ", ABS( S(iMaxError) - S_T ) / S(iMaxError)
  PRINT*, "CPU MINVAL(T) = ", MINVAL( T_S )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE FROM( T_S, Error_S )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE HOST( T_S, Error_S )
#endif
  DO iP = 1, nPoints
    IF( Error_S(iP) .NE. 0 ) CALL DescribeEOSInversionError( Error_S(iP) )
  END DO
  Error = ABS( T - T_S ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_S(iMaxError), Y(iMaxError), Ds_T, Ts_T, Ys_T, OS_S, Ss_T, S_T )
  PRINT*, "GPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "GPU  Error(S) = ", ABS( S(iMaxError) - S_T ) / S(iMaxError)
  PRINT*, "GPU MINVAL(T) = ", MINVAL( T_S )

  ! -------------------------------------------------------------------

  T_S = 0.0_dp
  Error_S = 0

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE TO( T_S, Error_S )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE DEVICE( T_S, Error_S )
#endif

  CALL CPU_TIME( tBegin )

  CALL ComputeTemperatureWith_DSY_Many &
         ( D, S, Y, Ds_T, Ts_T, Ys_T, Ss_T, OS_S, T_S, &
           UseInitialGuess_Option = .FALSE., &
           Error_Option = Error_S )

  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  CALL CPU_TIME( tBegin )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  !$OMP PRIVATE( T_Guess )
#elif defined (WEAKLIB_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR &
  !$ACC PRIVATE( T_Guess ) &
  !$ACC PRESENT( D, S, Y, Ds_T, Ts_T, Ys_T, Ss_T, OS_S, T_S, Error_S )
#elif defined (WEAKLIB_OMP)
  !$OMP PARALLEL DO &
  !$OMP PRIVATE( T_Guess )
#endif
  DO iP = 1, nPoints
    T_Guess = T_S(iP)
    CALL ComputeTemperatureWith_DSY_Single_NoGuess &
           ( D(iP), S(iP), Y(iP), Ds_T, Ts_T, Ys_T, Ss_T, OS_S, T_S(iP), &
             Error_Option = Error_S(iP) )
  END DO

  CALL CPU_TIME( tEnd )
  tGPU = tEnd - tBegin

  PRINT*
  PRINT*, "ComputeTemperatureWith_DSY (No Guess):"
  PRINT*, "CPU_TIME = ", tCPU, "GPU_TIME = ", tGPU

  DO iP = 1, nPoints
    IF( Error_S(iP) .NE. 0 ) CALL DescribeEOSInversionError( Error_S(iP) )
  END DO
  Error = ABS( T - T_S ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_S(iMaxError), Y(iMaxError), Ds_T, Ts_T, Ys_T, OS_S, Ss_T, S_T )
  PRINT*, "CPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "CPU  Error(S) = ", ABS( S(iMaxError) - S_T ) / S(iMaxError)
  PRINT*, "CPU MINVAL(T) = ", MINVAL( T_S )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE FROM( T_S, Error_S )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE HOST( T_S, Error_S )
#endif
  DO iP = 1, nPoints
    IF( Error_S(iP) .NE. 0 ) CALL DescribeEOSInversionError( Error_S(iP) )
  END DO
  Error = ABS( T - T_S ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_S(iMaxError), Y(iMaxError), Ds_T, Ts_T, Ys_T, OS_S, Ss_T, S_T )
  PRINT*, "GPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "GPU  Error(S) = ", ABS( S(iMaxError) - S_T ) / S(iMaxError)
  PRINT*, "GPU MINVAL(T) = ", MINVAL( T_S )

  ! -------------------------------------------------------------------

  ! -------------------------------------------------------------------
  ! --- Recover Temperature from Pressure ----------------------
  ! -------------------------------------------------------------------

  Amp = 0.001_dp
  dT = Amp * 2.0_dp * ( rndm_T - 0.5_dp ) * T

  associate &
    ( minT => 1.0001 * EOS % TS % MinValues(iT), &
      maxT => 0.9999 * EOS % TS % MaxValues(iT) )

  T_P = MIN( MAX( T + dT, minT ), maxT )

  end associate

  Error_P = 0

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE TO( T_P, Error_P )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE DEVICE( T_P, Error_P )
#endif

  CALL CPU_TIME( tBegin )

  CALL ComputeTemperatureWith_DPY_Many &
         ( D, P, Y, Ds_T, Ts_T, Ys_T, Ps_T, OS_P, T_P, &
           UseInitialGuess_Option = .TRUE., &
           Error_Option = Error_P )

  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  CALL CPU_TIME( tBegin )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  !$OMP PRIVATE( T_Guess )
#elif defined (WEAKLIB_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR &
  !$ACC PRIVATE( T_Guess ) &
  !$ACC PRESENT( D, P, Y, Ds_T, Ts_T, Ys_T, Ps_T, OS_P, T_P, Error_P )
#elif defined(WEAKLIB_OMP_OL)
  !$OMP PARALLEL DO &
  !$OMP PRIVATE( T_Guess )
#endif
  DO iP = 1, nPoints
    T_Guess = T_P(iP)
    CALL ComputeTemperatureWith_DPY_Single_Guess &
           ( D(iP), P(iP), Y(iP), Ds_T, Ts_T, Ys_T, Ps_T, OS_P, T_P(iP), T_Guess, &
             Error_Option = Error_P(iP) )
  END DO

  CALL CPU_TIME( tEnd )
  tGPU = tEnd - tBegin

  PRINT*
  PRINT*, "ComputeTemperatureWith_DPY (Good Guess):"
  PRINT*, "CPU_TIME = ", tCPU, "GPU_TIME = ", tGPU

  DO iP = 1, nPoints
    IF( Error_P(iP) .NE. 0 ) CALL DescribeEOSInversionError( Error_P(iP) )
  END DO
  Error = ABS( T - T_P ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_P(iMaxError), Y(iMaxError), Ds_T, Ts_T, Ys_T, OS_P, Ps_T, P_T )
  PRINT*, "CPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "CPU  Error(P) = ", ABS( P(iMaxError) - P_T ) / P(iMaxError)
  PRINT*, "CPU MINVAL(T) = ", MINVAL( T_P )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE FROM( T_P, Error_P )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE HOST( T_P, Error_P )
#endif
  DO iP = 1, nPoints
    IF( Error_P(iP) .NE. 0 ) CALL DescribeEOSInversionError( Error_P(iP) )
  END DO
  Error = ABS( T - T_P ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_P(iMaxError), Y(iMaxError), Ds_T, Ts_T, Ys_T, OS_P, Ps_T, P_T )
  PRINT*, "GPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "GPU  Error(P) = ", ABS( P(iMaxError) - P_T ) / P(iMaxError)
  PRINT*, "GPU MINVAL(T) = ", MINVAL( T_P )

  ! -------------------------------------------------------------------

  associate &
    ( minT => 1.0001 * EOS % TS % MinValues(iT), &
      maxT => 0.9999 * EOS % TS % MaxValues(iT) )

  Amp = 0.1_dp * ( maxT - minT ) ! --- Perturbation Amplitude
  dT = Amp * 2.0_dp * ( rndm_T - 0.5_dp ) * T

  T_P = MIN( MAX( T + dT, minT ), maxT )

  end associate

  Error_P = 0

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE TO( T_P, Error_P )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE DEVICE( T_P, Error_P )
#endif

  CALL CPU_TIME( tBegin )

  CALL ComputeTemperatureWith_DPY_Many &
         ( D, P, Y, Ds_T, Ts_T, Ys_T, Ps_T, OS_P, T_P, &
           UseInitialGuess_Option = .TRUE., &
           Error_Option = Error_P )

  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  CALL CPU_TIME( tBegin )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  !$OMP PRIVATE( T_Guess )
#elif defined (WEAKLIB_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR &
  !$ACC PRIVATE( T_Guess ) &
  !$ACC PRESENT( D, P, Y, Ds_T, Ts_T, Ys_T, Ps_T, OS_P, T_P, Error_P )
#elif defined(WEAKLIB_OMP_OL)
  !$OMP PARALLEL DO &
  !$OMP PRIVATE( T_Guess )
#endif
  DO iP = 1, nPoints
    T_Guess = T_P(iP)
    CALL ComputeTemperatureWith_DPY_Single_Guess &
           ( D(iP), P(iP), Y(iP), Ds_T, Ts_T, Ys_T, Ps_T, OS_P, T_P(iP), T_Guess, &
             Error_Option = Error_P(iP) )
  END DO

  CALL CPU_TIME( tEnd )
  tGPU = tEnd - tBegin

  PRINT*
  PRINT*, "ComputeTemperatureWith_DPY (Bad Guess):"
  PRINT*, "CPU_TIME = ", tCPU, "GPU_TIME = ", tGPU

  DO iP = 1, nPoints
    IF( Error_P(iP) .NE. 0 ) CALL DescribeEOSInversionError( Error_P(iP) )
  END DO
  Error = ABS( T - T_P ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_P(iMaxError), Y(iMaxError), Ds_T, Ts_T, Ys_T, OS_P, Ps_T, P_T )
  PRINT*, "CPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "CPU  Error(P) = ", ABS( P(iMaxError) - P_T ) / P(iMaxError)
  PRINT*, "CPU MINVAL(T) = ", MINVAL( T_P )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE FROM( T_P, Error_P )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE HOST( T_P, Error_P )
#endif
  DO iP = 1, nPoints
    IF( Error_P(iP) .NE. 0 ) CALL DescribeEOSInversionError( Error_P(iP) )
  END DO
  Error = ABS( T - T_P ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_P(iMaxError), Y(iMaxError), Ds_T, Ts_T, Ys_T, OS_P, Ps_T, P_T )
  PRINT*, "GPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "GPU  Error(P) = ", ABS( P(iMaxError) - P_T ) / P(iMaxError)
  PRINT*, "GPU MINVAL(T) = ", MINVAL( T_P )

  ! -------------------------------------------------------------------

  T_P = 0.0_dp
  Error_P = 0

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE TO( T_P, Error_P )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE DEVICE( T_P, Error_P )
#endif

  CALL CPU_TIME( tBegin )

  CALL ComputeTemperatureWith_DPY_Many &
         ( D, P, Y, Ds_T, Ts_T, Ys_T, Ps_T, OS_P, T_P, &
           UseInitialGuess_Option = .FALSE., &
           Error_Option = Error_P )

  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  CALL CPU_TIME( tBegin )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  !$OMP PRIVATE( T_Guess )
#elif defined (WEAKLIB_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR &
  !$ACC PRIVATE( T_Guess ) &
  !$ACC PRESENT( D, P, Y, Ds_T, Ts_T, Ys_T, Ps_T, OS_P, T_P, Error_P )
#elif defined (WEAKLIB_OMP)
  !$OMP PARALLEL DO &
  !$OMP PRIVATE( T_Guess )
#endif
  DO iP = 1, nPoints
    T_Guess = T_P(iP)
    CALL ComputeTemperatureWith_DPY_Single_NoGuess &
           ( D(iP), P(iP), Y(iP), Ds_T, Ts_T, Ys_T, Ps_T, OS_P, T_P(iP), &
             Error_Option = Error_P(iP) )
  END DO

  CALL CPU_TIME( tEnd )
  tGPU = tEnd - tBegin

  PRINT*
  PRINT*, "ComputeTemperatureWith_DPY (No Guess):"
  PRINT*, "CPU_TIME = ", tCPU, "GPU_TIME = ", tGPU

  DO iP = 1, nPoints
    IF( Error_P(iP) .NE. 0 ) CALL DescribeEOSInversionError( Error_P(iP) )
  END DO
  Error = ABS( T - T_P ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_P(iMaxError), Y(iMaxError), Ds_T, Ts_T, Ys_T, OS_P, Ps_T, P_T )
  PRINT*, "CPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "CPU  Error(P) = ", ABS( P(iMaxError) - P_T ) / P(iMaxError)
  PRINT*, "CPU MINVAL(T) = ", MINVAL( T_P )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE FROM( T_P, Error_P )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE HOST( T_P, Error_P )
#endif
  DO iP = 1, nPoints
    IF( Error_P(iP) .NE. 0 ) CALL DescribeEOSInversionError( Error_P(iP) )
  END DO
  Error = ABS( T - T_P ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_P(iMaxError), Y(iMaxError), Ds_T, Ts_T, Ys_T, OS_P, Ps_T, P_T )
  PRINT*, "GPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "GPU  Error(P) = ", ABS( P(iMaxError) - P_T ) / P(iMaxError)
  PRINT*, "GPU MINVAL(T) = ", MINVAL( T_P )

  ! -------------------------------------------------------------------

  !CALL CPU_TIME( tBegin )

  !DO iP = 1, nPoints

  !  CALL ComputeTempFromEntropy &
  !         ( D(iP), S(iP), Y(iP), Ds_T, Ts_T, Ys_T, [1,1,0], &
  !           Ss_T, OS_S, T_S(iP:iP) )

  !END DO

  !CALL CPU_TIME( tEnd )

  !Error = ABS( T - T_S ) / T

  !PRINT*
  !PRINT*, "ComputeTempFromEntropy:"
  !PRINT*, "CPU_TIME = ", tEnd - tBegin
  !PRINT*, MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  !PRINT*, "MINVAL(T) = ", MINVAL( T_S )

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
  !!$ACC SHUTDOWN
#endif

  CALL MPI_FINALIZE( ierr )

END PROGRAM wlEosInversionTest
