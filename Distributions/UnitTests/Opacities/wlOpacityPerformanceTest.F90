PROGRAM wlOpacityPerformanceTest

  USE wlKindModule, ONLY: &
    dp
  USE wlOpacityTableModule, ONLY: &
    OpacityTableType
  USE wlIOModuleHDF, ONLY: &
    InitializeHDF, &
    FinalizeHDF
  USE wlOpacityTableIOModuleHDF, ONLY: &
    ReadOpacityTableHDF
  USE wlInterpolationModule, ONLY: &
    LogInterpolateSingleVariable_1D3D_Custom, &
    LogInterpolateSingleVariable_1D3D_Custom_Point, &
    LogInterpolateSingleVariable_4D_Custom_Point

  IMPLICIT NONE

  INTEGER :: &
    n_rndm, &
    iP, jP
  INTEGER, PARAMETER :: &
    iD = 1, iT = 2, iY = 3, &
    nPointsX = 2**18, &
    nPointsE = 2**05
  REAL(dp) :: &
    tBegin, &
    tEnd, &
    tCPU, &
    tGPU
  REAL(dp), DIMENSION(:), ALLOCATABLE :: &
    LogEs_T, LogDs_T, LogTs_T
  REAL(dp), DIMENSION(nPointsX) :: &
    Y, &
    D, LogD, &
    T, LogT, &
    rndm_D, &
    rndm_T, &
    rndm_Y
  REAL(dp), DIMENSION(nPointsE) :: &
    E, LogE, &
    rndm_E
  REAL(dp), DIMENSION(nPointsE,nPointsX) :: &
    opEC_CPU_1, opEC_GPU_1, &
    opEC_CPU_2, opEC_GPU_2, &
    opEC_CPU_3, opEC_GPU_3
  TYPE(OpacityTableType) :: &
    OpTab

#if defined(WEAKLIB_OMP_OL)
#elif defined(WEAKLIB_OACC)
  !$ACC INIT
#endif

  CALL InitializeHDF( )

  CALL ReadOpacityTableHDF &
         ( OpTab, "wl-Op-SFHo-15-25-50-E40-B85-AbEm.h5" )

  CALL FinalizeHDF( )

  ASSOCIATE &
    ( Es_T    => OpTab % EnergyGrid % Values,                  &
      Ds_T    => OpTab % EOSTable % TS % States(iD)  % Values, &
      Ts_T    => OpTab % EOSTable % TS % States(iT)  % Values, &
      Ys_T    => OpTab % EOSTable % TS % States(iY)  % Values, &
      EmAb_T  => OpTab % EmAb % Opacity(1) % Values, &
      OS_EmAb => OpTab % EmAb % Offsets(1) )

  ALLOCATE( LogEs_T(SIZE( Es_T )) )
  LogEs_T = LOG10( Es_T )

  ALLOCATE( LogDs_T(SIZE( Ds_T )) )
  LogDs_T = LOG10( Ds_T )

  ALLOCATE( LogTs_T(SIZE( Ts_T )) )
  LogTs_T = LOG10( Ts_T )

  WRITE(*,*)
  WRITE(*,'(A4,A12,ES10.4E2,A3,ES10.4E2)') &
    '', 'Max/Min E = ', OpTab % EnergyGrid % MinValue, &
    ' / ', OpTab % EnergyGrid % MaxValue
  WRITE(*,'(A4,A12,ES10.4E2,A3,ES10.4E2)') &
    '', 'Max/Min D = ', OpTab % EOSTable % TS % MinValues(iD), &
    ' / ', OpTab % EOSTable % TS % MaxValues(iD)
  WRITE(*,'(A4,A12,ES10.4E2,A3,ES10.4E2)') &
    '', 'Max/Min T = ', OpTab % EOSTable % TS % MinValues(iT), &
    ' / ', OpTab % EOSTable % TS % MaxValues(iT)
  WRITE(*,'(A4,A12,ES10.4E2,A3,ES10.4E2)') &
    '', 'Max/Min Y = ', OpTab % EOSTable % TS % MinValues(iY), &
    ' / ', OpTab % EOSTable % TS % MaxValues(iY)
  WRITE(*,*)
  WRITE(*,'(A4,A14,I10.10,A7)') &
    '', 'Interpolating ', nPointsE * nPointsX, ' points'
  WRITE(*,*)

  n_rndm = nPointsE

  ! --- Initialize Energy Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_E )

  ASSOCIATE &
    ( minE => OpTab % EnergyGrid % MinValue, &
      maxE => OpTab % EnergyGrid % MaxValue )

  E(:) = 10**( LOG10(minE) + ( LOG10(maxE) - LOG10(minE) ) * rndm_E )

  END ASSOCIATE

  n_rndm = nPointsX

  ! --- Initialize Density Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_D )

  ASSOCIATE &
    ( minD => OpTab % EOSTable % TS % MinValues(iD), &
      maxD => OpTab % EOSTable % TS % MaxValues(iD) )

  D(:) = 10**( LOG10(minD) + ( LOG10(maxD) - LOG10(minD) ) * rndm_D )

  END ASSOCIATE

  ! --- Initialize Temperature Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_T )

  ASSOCIATE &
    ( minT => OpTab % EOSTable % TS % MinValues(iT), &
      maxT => OpTab % EOSTable % TS % MaxValues(iT) )

  T(:) = 10**( LOG10(minT) + ( LOG10(maxT) - LOG10(minT) ) * rndm_T )

  END ASSOCIATE

  ! --- Initialize Electron Fraction Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_Y )

  ASSOCIATE &
    ( minY => OpTab % EOSTable % TS % MinValues(iY), &
      maxY => OpTab % EOSTable % TS % MaxValues(iY) )

  Y(:) = minY + ( maxY - minY ) * rndm_Y

  END ASSOCIATE

  LogE = LOG10( E )
  LogD = LOG10( D )
  LogT = LOG10( T )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET ENTER DATA &
  !$OMP MAP( to: LogEs_T, LogDs_T, LogTs_T, Ys_T, OS_EmAb, LogE, LogD, LogT, Y, EmAb_T ) &
  !$OMP MAP( alloc: opEC_GPU_1, opEC_GPU_2, opEC_GPU_3 )
#elif defined (WEAKLIB_OACC)
  !$ACC ENTER DATA &
  !$ACC COPYIN( LogEs_T, LogDs_T, LogTs_T, Ys_T, OS_EmAb, LogE, LogD, LogT, Y, EmAb_T ) &
  !$ACC CREATE( opEC_GPU_1, opEC_GPU_2, opEC_GPU_3 )
#endif

  tBegin = get_wtime()
  CALL LogInterpolateSingleVariable_1D3D_Custom &
         ( LogE, LogD, LogT, Y, LogEs_T, LogDs_T, LogTs_T, Ys_T, OS_EmAb, EmAb_T, &
           opEC_CPU_1, GPU_Option = .FALSE. )
  tEnd = get_wtime()
  tCPU = tEnd - tBegin

  tBegin = get_wtime()
  CALL LogInterpolateSingleVariable_1D3D_Custom &
         ( LogE, LogD, LogT, Y, LogEs_T, LogDs_T, LogTs_T, Ys_T, OS_EmAb, EmAb_T, &
           opEC_GPU_1, GPU_Option = .TRUE. )
  tEnd = get_wtime()
  tGPU = tEnd - tBegin

  WRITE(*,*)
  WRITE(*,'(A4,A60,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_1D3D_Custom (CPU): ', tCPU
  WRITE(*,'(A4,A60,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_1D3D_Custom (GPU): ', tGPU
  WRITE(*,*)

  tBegin = get_wtime()
  DO iP = 1, nPointsX

    CALL LogInterpolateSingleVariable_1D3D_Custom_Point &
           ( LogE,    LogD(iP), LogT(iP), Y(iP), &
             LogEs_T, LogDs_T,  LogTs_T,  Ys_T,  &
             OS_EmAb, EmAb_T, opEC_CPU_2(:,iP) )

  END DO
  tEnd = get_wtime()
  tCPU = tEnd - tBegin

  tBegin = get_wtime()
#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined (WEAKLIB_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR &
  !$ACC PRESENT( LogE, LogD, LogT, Y, LogEs_T, LogDs_T, LogTs_T, Ys_T, &
  !$ACC          OS_EmAb, EmAb_T, opEC_GPU_2 )
#endif
  DO iP = 1, nPointsX

    CALL LogInterpolateSingleVariable_1D3D_Custom_Point &
           ( LogE,    LogD(iP), LogT(iP), Y(iP), &
             LogEs_T, LogDs_T,  LogTs_T,  Ys_T,  &
             OS_EmAb, EmAb_T, opEC_GPU_2(:,iP) )

  END DO
  tEnd = get_wtime()
  tGPU = tEnd - tBegin

  WRITE(*,*)
  WRITE(*,'(A4,A60,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_1D3D_Custom_Point (CPU): ', tCPU
  WRITE(*,'(A4,A60,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_1D3D_Custom_Point (GPU): ', tGPU
  WRITE(*,*)

  tBegin = get_wtime()
  DO jP = 1, nPointsX
    DO iP = 1, nPointsE

      CALL LogInterpolateSingleVariable_4D_Custom_Point &
             ( LogE(iP), LogD(jP), LogT(jP), Y(jP), &
               LogEs_T,  LogDs_T,  LogTs_T,  Ys_T,  &
               OS_EmAb, EmAb_T, opEC_CPU_3(iP,jP) )

    END DO
  END DO
  tEnd = get_wtime()
  tCPU = tEnd - tBegin

  tBegin = get_wtime()
#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(2)
#elif defined (WEAKLIB_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(2) &
  !$ACC PRESENT( LogE, LogD, LogT, Y, LogEs_T, LogDs_T, LogTs_T, Ys_T, &
  !$ACC          OS_EmAb, EmAb_T, opEC_GPU_3 )
#endif
  DO jP = 1, nPointsX
    DO iP = 1, nPointsE

      CALL LogInterpolateSingleVariable_4D_Custom_Point &
             ( LogE(iP), LogD(jP), LogT(jP), Y(jP), &
               LogEs_T,  LogDs_T,  LogTs_T,  Ys_T,  &
               OS_EmAb, EmAb_T, opEC_GPU_3(iP,jP) )

    END DO
  END DO
  tEnd = get_wtime()
  tGPU = tEnd - tBegin

  WRITE(*,*)
  WRITE(*,'(A4,A60,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_4D_Custom_Point (CPU): ', tCPU
  WRITE(*,'(A4,A60,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_4D_Custom_Point (GPU): ', tGPU
  WRITE(*,*)

  WRITE(*,*) 'CPU Results'
  WRITE(*,*) MINVAL( opEC_CPU_1 ), MAXVAL( opEC_CPU_1 )
  WRITE(*,*) MINVAL( opEC_CPU_2 ), MAXVAL( opEC_CPU_2 )
  WRITE(*,*) MINVAL( opEC_CPU_3 ), MAXVAL( opEC_CPU_3 )
  WRITE(*,*) MAXVAL( ABS( opEC_CPU_1 - opEC_CPU_2 ) )
  WRITE(*,*) MAXVAL( ABS( opEC_CPU_1 - opEC_CPU_3 ) )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE FROM( opEC_GPU_1, opEC_GPU_2, opEC_GPU_3 )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE HOST( opEC_GPU_1, opEC_GPU_2, opEC_GPU_3 )
#endif

  WRITE(*,*) 'GPU Results'
  WRITE(*,*) MINVAL( opEC_GPU_1 ), MAXVAL( opEC_GPU_1 )
  WRITE(*,*) MINVAL( opEC_GPU_2 ), MAXVAL( opEC_GPU_2 )
  WRITE(*,*) MINVAL( opEC_GPU_3 ), MAXVAL( opEC_GPU_3 )
  WRITE(*,*) MAXVAL( ABS( opEC_GPU_1 - opEC_GPU_2 ) )
  WRITE(*,*) MAXVAL( ABS( opEC_GPU_1 - opEC_GPU_3 ) )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET EXIT DATA &
  !$OMP MAP( release: LogEs_T, LogDs_T, LogTs_T, Ys_T, OS_EmAb, LogE, LogD, LogT, Y, EmAb_T, &
  !$OMP               opEC_GPU_1, opEC_GPU_2, opEC_GPU_3 )
#elif defined (WEAKLIB_OACC)
  !$ACC EXIT DATA &
  !$ACC DELETE( LogEs_T, LogDs_T, LogTs_T, Ys_T, OS_EmAb, LogE, LogD, LogT, Y, EmAb_T )
#endif

  END ASSOCIATE

#if defined(WEAKLIB_OMP_OL)
#elif defined(WEAKLIB_OACC)
  !$ACC SHUTDOWN
#endif

CONTAINS

  REAL(dp) FUNCTION get_wtime()
    USE, INTRINSIC :: iso_fortran_env, only: i8=>int64, dp=>real64
    IMPLICIT NONE

    INTEGER(i8) :: clock_read
    INTEGER(i8) :: clock_rate
    INTEGER(i8) :: clock_max

    CALL SYSTEM_CLOCK(clock_read,clock_rate,clock_max)
    get_wtime = REAL(clock_read,dp) / REAL(clock_rate,dp)

    RETURN
  END FUNCTION get_wtime


END PROGRAM wlOpacityPerformanceTest
