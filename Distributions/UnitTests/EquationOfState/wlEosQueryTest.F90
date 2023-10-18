PROGRAM wlEosQueryTest

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

  IMPLICIT NONE

  REAL(dp), PARAMETER :: D_0 = 1660.55
  REAL(dp), PARAMETER :: T_0 = 1.16046e9
  REAL(dp), PARAMETER :: Y_0 = 0.5

  INTEGER :: &
    iP
  INTEGER, PARAMETER :: &
    nPoints = 2**16, &
    iD = 1, iT = 2, iY = 3
  REAL(dp) :: &
    tBegin, &
    tEnd, &
    tCPU, &
    tGPU
  REAL(dp), DIMENSION(nPoints) :: &
    D, T, Y, E, P, S, ET
  TYPE(EquationOfStateTableType) :: &
    EOS
  REAL(DP), DIMENSION(:), ALLOCATABLE :: &
    Ds_T, Ts_T, Ys_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: &
    Ps_T, Es_T, Ss_T, ETs_T
  REAL(DP) :: &
    OS_P, OS_E, OS_S, OS_ET
  INTEGER :: &
    iD_T, iT_T, iY_T, iP_T, iS_T, iE_T, iET_T

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

  iP_T  = EOS % DV % Indices % iPressure
  iS_T  = EOS % DV % Indices % iEntropyPerBaryon
  iE_T  = EOS % DV % Indices % iInternalEnergyDensity
  iET_T = EOS % DV % Indices % iThermalEnergy

  OS_P  =  EOS % DV % Offsets(iP_T)
  OS_S  = EOS % DV % Offsets(iS_T)
  OS_E  = EOS % DV % Offsets(iE_T)
  OS_ET = EOS % DV % Offsets(iET_T)

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
  ALLOCATE &
    ( ETs_T(1:EOS % DV % nPoints(1), &
            1:EOS % DV % nPoints(2), &
            1:EOS % DV % nPoints(3)) )

  Ps_T  = EOS % DV % Variables(iP_T )  % Values
  Ss_T  = EOS % DV % Variables(iS_T )  % Values
  Es_T  = EOS % DV % Variables(iE_T )  % Values
  ETs_T = EOS % DV % Variables(iET_T ) % Values

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET ENTER DATA &
  !$OMP MAP( to: Ds_T, Ts_T, Ys_T, Es_T, ETs_T, Ps_T, Ss_T, OS_E, OS_P, OS_S, OS_ET ) &
  !$OMP MAP( alloc: D, T, Y, P, E, S, ET )
#elif defined (WEAKLIB_OACC)
  !$ACC ENTER DATA &
  !$ACC COPYIN( Ds_T, Ts_T, Ys_T, Es_T, ETs_T, Ps_T, Ss_T, OS_E, OS_P, OS_S, OS_ET ) &
  !$ACC CREATE( D, T, Y, P, E, S, ET )
#endif

  WRITE(*,*)
  WRITE(*,'(A4,A10,I10.10)') '', 'nPoints = ', nPoints

  ! --- Initialize Density Points ---

  D = D_0


  ! --- Initialize Temperature Points ---

  T = T_0


  ! --- Initialize Electron Fraction Points ---

  Y = Y_0

  PRINT*, "Min/Max D = ", MINVAL( D ), MAXVAL( D )
  PRINT*, "Min/Max T = ", MINVAL( T ), MAXVAL( T )
  PRINT*, "Min/Max Y = ", MINVAL( Y ), MAXVAL( Y )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE TO( D, Y, T )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE DEVICE( D, Y, T )
#endif

  ! --- Compute Internal Energy, Pressure, and Entropy ---

  DO iP = 1, nPoints

    CALL LogInterpolateSingleVariable_3D_Custom_Point &
           ( D(iP), T(iP), Y(iP), Ds_T, Ts_T, Ys_T, OS_E, Es_T, E(iP) )

    CALL LogInterpolateSingleVariable_3D_Custom_Point &
           ( D(iP), T(iP), Y(iP), Ds_T, Ts_T, Ys_T, OS_P, Ps_T, P(iP) )

    CALL LogInterpolateSingleVariable_3D_Custom_Point &
           ( D(iP), T(iP), Y(iP), Ds_T, Ts_T, Ys_T, OS_S, Ss_T, S(iP) )

    CALL LogInterpolateSingleVariable_3D_Custom_Point &
           ( D(iP), T(iP), Y(iP), Ds_T, Ts_T, Ys_T, OS_ET, ETs_T, ET(iP) )

  END DO

  WRITE(*,*)
  WRITE(*,*) "Min/Max E = ", MINVAL( E ), MAXVAL( E )
  WRITE(*,*) "Min/Max P = ", MINVAL( P ), MAXVAL( P )
  WRITE(*,*) "Min/Max S = ", MINVAL( S ), MAXVAL( S )
  WRITE(*,*) "Min/Max ET = ", MINVAL( ET ), MAXVAL( ET )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET EXIT DATA &
  !$OMP MAP( release: D, E, ET, P, S, Y, Ds_T, Ts_T, Ys_T, Es_T, ETs_T, Ps_T, Ss_T, OS_E, OS_P, OS_S, OS_ET )
#elif defined (WEAKLIB_OACC)
  !$ACC EXIT DATA &
  !$ACC DELETE( D, E, ET, P, S, Y, Ds_T, Ts_T, Ys_T, Es_T, ETs_T, Ps_T, Ss_T, OS_E, OS_P, OS_S, OS_ET )
#endif

  DEALLOCATE( Ds_T, Ts_T, Ys_T, Ps_T, Ss_T, Es_T )

#if defined(WEAKLIB_OMP_OL)
#elif defined(WEAKLIB_OACC)
  !$ACC SHUTDOWN
#endif

END PROGRAM wlEosQueryTest
