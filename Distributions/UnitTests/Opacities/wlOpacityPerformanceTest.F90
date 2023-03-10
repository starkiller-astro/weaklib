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
    LogInterpolateSingleVariable_2D_Custom_Point, &
    LogInterpolateSingleVariable_3D_Custom, &
    LogInterpolateSingleVariable_1D3D_Custom, &
    LogInterpolateSingleVariable_1D3D_Custom_Point, &
    LogInterpolateSingleVariable_2D2D_Custom, &
    LogInterpolateSingleVariable_2D2D_Custom_Point, &
    LogInterpolateSingleVariable_2D2D_Custom_Aligned, &
    LogInterpolateSingleVariable_2D2D_Custom_Aligned_Point, &
    LogInterpolateSingleVariable_4D_Custom_Point

  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: &
    n_rndm, &
    iP, jP, &
    iM, iEta, iT, iE2, iE1, &
    nMoments, nPointsEta, nPointsT, &
    iH1, iH2, &
    ierr
  INTEGER, PARAMETER :: &
    iD_T = 1, iT_T = 2, iY_T = 3, &
    nPointsX = 2**15, &
    nPointsE = 2**05
  REAL(dp) :: &
    tBegin, &
    tEnd, &
    tCPU, &
    tGPU
  REAL(dp), PARAMETER :: &
    kMeV = 8.61733e-11_dp
  REAL(dp), DIMENSION(nPointsX) :: &
    Y, &
    D, LogD, &
    T, LogT, &
    Eta, LogEta, &
    rndm_D, &
    rndm_T, &
    rndm_Y, &
    rndm_Eta
  REAL(dp), DIMENSION(nPointsE) :: &
    E, LogE, &
    rndm_E
  REAL(dp), DIMENSION(:), ALLOCATABLE :: &
    LogEs_T, LogDs_T, LogTs_T, LogEtas_T
  REAL(dp), DIMENSION(:,:,:,:,:), ALLOCATABLE :: &
    NES_T, NES_AT
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: &
    opEC_CPU_1, opEC_GPU_1, &
    opEC_CPU_2, opEC_GPU_2, &
    opEC_CPU_3, opEC_GPU_3
  REAL(dp), DIMENSION(:,:,:), ALLOCATABLE :: &
    opH1_CPU_1, opH1_GPU_1, &
    opH1_CPU_2, opH1_GPU_2, &
    opH1_CPU_3, opH1_GPU_3, &
    opH1_CPU_4, opH1_GPU_4
  TYPE(OpacityTableType) :: &
    OpTab

  CALL MPI_INIT( ierr )

#if defined(WEAKLIB_OMP_OL)
#elif defined(WEAKLIB_OACC)
  !!$ACC INIT
#endif

  CALL InitializeHDF( )

  CALL ReadOpacityTableHDF &
         ( OpTab, &
           FileName_EmAb_Option = "wl-Op-SFHo-15-25-50-E40-B85-AbEm.h5", &
           FileName_Iso_Option  = "wl-Op-SFHo-15-25-50-E40-B85-Iso.h5", &
           FileName_NES_Option  = "wl-Op-SFHo-15-25-50-E40-B85-NES.h5", &
           FileName_Pair_Option = "wl-Op-SFHo-15-25-50-E40-B85-Pair.h5", &
           EquationOfStateTableName_Option = "EquationOfStateTable.h5" )

  CALL FinalizeHDF( )

  ALLOCATE( opEC_CPU_1(nPointsE,nPointsX) )
  ALLOCATE( opEC_CPU_2(nPointsE,nPointsX) )
  ALLOCATE( opEC_CPU_3(nPointsE,nPointsX) )
  ALLOCATE( opEC_GPU_1(nPointsE,nPointsX) )
  ALLOCATE( opEC_GPU_2(nPointsE,nPointsX) )
  ALLOCATE( opEC_GPU_3(nPointsE,nPointsX) )

  ALLOCATE( opH1_CPU_1(nPointsE,nPointsE,nPointsX) )
  ALLOCATE( opH1_CPU_2(nPointsE,nPointsE,nPointsX) )
  ALLOCATE( opH1_CPU_3(nPointsE,nPointsE,nPointsX) )
  ALLOCATE( opH1_CPU_4(nPointsE,nPointsE,nPointsX) )
  ALLOCATE( opH1_GPU_1(nPointsE,nPointsE,nPointsX) )
  ALLOCATE( opH1_GPU_2(nPointsE,nPointsE,nPointsX) )
  ALLOCATE( opH1_GPU_3(nPointsE,nPointsE,nPointsX) )
  ALLOCATE( opH1_GPU_4(nPointsE,nPointsE,nPointsX) )

  ASSOCIATE &
    ( Es_T    => OpTab % EnergyGrid % Values, &
      Ds_T    => OpTab % EOSTable % TS % States(iD_T) % Values, &
      Ts_T    => OpTab % EOSTable % TS % States(iT_T) % Values, &
      Ys_T    => OpTab % EOSTable % TS % States(iY_T) % Values, &
      Etas_T  => OpTab % EtaGrid % Values, &
      EmAb_T  => OpTab % EmAb % Opacity(1) % Values, &
      OS_EmAb => OpTab % EmAb % Offsets(1), &
      OS_NES  => OpTab % Scat_NES % Offsets )

  ALLOCATE( LogEs_T(SIZE( Es_T )) )
  LogEs_T = LOG10( Es_T )

  ALLOCATE( LogDs_T(SIZE( Ds_T )) )
  LogDs_T = LOG10( Ds_T )

  ALLOCATE( LogTs_T(SIZE( Ts_T )) )
  LogTs_T = LOG10( Ts_T )

  ALLOCATE( LogEtas_T(SIZE( Etas_T )) )
  LogEtas_T = LOG10( Etas_T )

  ALLOCATE( NES_T(OpTab % Scat_NES % nPoints(1), &
                  OpTab % Scat_NES % nPoints(2), &
                  OpTab % Scat_NES % nPoints(4), &
                  OpTab % Scat_NES % nPoints(5), &
                  OpTab % Scat_NES % nMoments) )
  DO iM = 1, OpTab % Scat_NES % nMoments
    NES_T(:,:,:,:,iM) = OpTab % Scat_NES % Kernel(1) % Values(:,:,iM,:,:)
  END DO

  ALLOCATE( NES_AT(nPointsE, &
                   nPointsE, &
                   OpTab % Scat_NES % nPoints(4), &
                   OpTab % Scat_NES % nPoints(5), &
                   OpTab % Scat_NES % nMoments) )
  NES_AT = 0.0d0

  WRITE(*,*)
  WRITE(*,'(A4,A14,ES10.4E2,A3,ES10.4E2)') &
    '', 'Max/Min E   = ', OpTab % EnergyGrid % MinValue, &
    ' / ', OpTab % EnergyGrid % MaxValue
  WRITE(*,'(A4,A14,ES10.4E2,A3,ES10.4E2)') &
    '', 'Max/Min D   = ', OpTab % EOSTable % TS % MinValues(iD_T), &
    ' / ', OpTab % EOSTable % TS % MaxValues(iD_T)
  WRITE(*,'(A4,A14,ES10.4E2,A3,ES10.4E2)') &
    '', 'Max/Min T   = ', OpTab % EOSTable % TS % MinValues(iT_T), &
    ' / ', OpTab % EOSTable % TS % MaxValues(iT_T)
  WRITE(*,'(A4,A14,ES10.4E2,A3,ES10.4E2)') &
    '', 'Max/Min Y   = ', OpTab % EOSTable % TS % MinValues(iY_T), &
    ' / ', OpTab % EOSTable % TS % MaxValues(iY_T)
  WRITE(*,'(A4,A14,ES10.4E2,A3,ES10.4E2)') &
    '', 'Max/Min Eta = ', OpTab % EtaGrid % MinValue, &
    ' / ', OpTab % EtaGrid % MaxValue
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
    ( minD => OpTab % EOSTable % TS % MinValues(iD_T), &
      maxD => OpTab % EOSTable % TS % MaxValues(iD_T) )

  D(:) = 10**( LOG10(minD) + ( LOG10(maxD) - LOG10(minD) ) * rndm_D )

  END ASSOCIATE

  ! --- Initialize Temperature Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_T )

  ASSOCIATE &
    ( minT => OpTab % EOSTable % TS % MinValues(iT_T), &
      maxT => OpTab % EOSTable % TS % MaxValues(iT_T) )

  T(:) = 10**( LOG10(minT) + ( LOG10(maxT) - LOG10(minT) ) * rndm_T )

  END ASSOCIATE

  ! --- Initialize Electron Fraction Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_Y )

  ASSOCIATE &
    ( minY => OpTab % EOSTable % TS % MinValues(iY_T), &
      maxY => OpTab % EOSTable % TS % MaxValues(iY_T) )

  Y(:) = minY + ( maxY - minY ) * rndm_Y

  END ASSOCIATE

  ! --- Initialize Eta Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_Eta )

  ASSOCIATE &
    ( minEta => OpTab % EtaGrid % MinValue, &
      maxEta => OpTab % EtaGrid % MaxValue )

  Eta(:) = minEta + ( maxEta - minEta ) * rndm_Eta

  END ASSOCIATE

  LogE = LOG10( E )
  LogD = LOG10( D )
  LogT = LOG10( T )
  LogEta = LOG10( Eta )

  opEC_GPU_1 = 0.0d0
  opEC_GPU_2 = 0.0d0
  opEC_GPU_3 = 0.0d0

  opH1_GPU_1 = 0.0d0
  opH1_GPU_2 = 0.0d0
  opH1_GPU_3 = 0.0d0
  opH1_GPU_4 = 0.0d0

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET ENTER DATA &
  !$OMP MAP( to: LogEs_T, LogDs_T, LogTs_T, Ys_T, LogEtas_T, &
  !$OMP          LogE, LogD, LogT, Y, LogEta, &
  !$OMP          OS_EmAb, EmAb_T, OS_NES, NES_T, NES_AT, &
  !$OMP          opEC_GPU_1, opEC_GPU_2, opEC_GPU_3, &
  !$OMP          opH1_GPU_1, opH1_GPU_2, opH1_GPU_3, opH1_GPU_4 )
#elif defined (WEAKLIB_OACC)
  !$ACC ENTER DATA &
  !$ACC COPYIN( LogEs_T, LogDs_T, LogTs_T, Ys_T, LogEtas_T, &
  !$ACC         LogE, LogD, LogT, Y, LogEta, &
  !$ACC         OS_EmAb, EmAb_T, OS_NES, NES_T, NES_AT, &
  !$ACC         opEC_GPU_1, opEC_GPU_2, opEC_GPU_3, &
  !$ACC         opH1_GPU_1, opH1_GPU_2, opH1_GPU_3, opH1_GPU_4 )
#endif

  nMoments   = OpTab % Scat_NES % nMoments
  nPointsEta = OpTab % Scat_NES % nPoints(5)
  nPointsT   = OpTab % Scat_NES % nPoints(4)

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5)
#elif defined (WEAKLIB_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
  !$ACC PRESENT( LogE, LogEs_T, OS_NES, NES_T, NES_AT )
#endif
  DO iM = 1, nMoments
    DO iEta = 1, nPointsEta
      DO iT = 1, nPointsT
        DO iE2 = 1, nPointsE
          DO iE1 = 1, nPointsE

            CALL LogInterpolateSingleVariable_2D_Custom_Point &
                   ( LogE(iE1), LogE(iE2), LogEs_T, LogEs_T, OS_NES(1,iM), &
                     NES_T(:,:,iT,iEta,iM), NES_AT(iE1,iE2,iT,iEta,iM) )

            NES_AT(iE1,iE2,iT,iEta,iM) &
              = LOG10( NES_AT(iE1,iE2,iT,iEta,iM) + OS_NES(1,iM) )

          END DO
        END DO
      END DO
    END DO
  END DO
#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE FROM( NES_AT )
#elif defined(WEAKLIB_OACC)
  !$ACC UPDATE HOST( NES_AT )
#endif

  !tBegin = get_wtime()
  !CALL LogInterpolateSingleVariable_1D3D_Custom &
  !       ( LogE, LogD, LogT, Y, LogEs_T, LogDs_T, LogTs_T, Ys_T, OS_EmAb, EmAb_T, &
  !         opEC_CPU_1 )
  !tEnd = get_wtime()
  tCPU = 0.0d0

  tBegin = get_wtime()
  CALL LogInterpolateSingleVariable_1D3D_Custom &
         ( LogE, LogD, LogT, Y, LogEs_T, LogDs_T, LogTs_T, Ys_T, OS_EmAb, EmAb_T, &
           opEC_GPU_1 )
  tEnd = get_wtime()
  tGPU = tEnd - tBegin

  WRITE(*,*)
  WRITE(*,'(A4,A64,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_1D3D_Custom (CPU): ', tCPU
  WRITE(*,'(A4,A64,ES10.4E2)') &
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
  opEC_CPU_1 = opEC_CPU_2

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
  WRITE(*,'(A4,A64,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_1D3D_Custom_Point (CPU): ', tCPU
  WRITE(*,'(A4,A64,ES10.4E2)') &
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
  WRITE(*,'(A4,A64,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_4D_Custom_Point (CPU): ', tCPU
  WRITE(*,'(A4,A64,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_4D_Custom_Point (GPU): ', tGPU
  WRITE(*,*)

  WRITE(*,*) 'EmAb (CPU)'
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

  WRITE(*,*) 'EmAb (GPU)'
  WRITE(*,*) MINVAL( opEC_GPU_1 ), MAXVAL( opEC_GPU_1 )
  WRITE(*,*) MINVAL( opEC_GPU_2 ), MAXVAL( opEC_GPU_2 )
  WRITE(*,*) MINVAL( opEC_GPU_3 ), MAXVAL( opEC_GPU_3 )
  WRITE(*,*) MAXVAL( ABS( opEC_GPU_1 - opEC_GPU_2 ) )
  WRITE(*,*) MAXVAL( ABS( opEC_GPU_1 - opEC_GPU_3 ) )

  iM = 1
  iH1 = ( iM - 1 ) * 2 + 1
  iH2 = ( iM - 1 ) * 2 + 2

  !tBegin = get_wtime()
  !CALL LogInterpolateSingleVariable_2D2D_Custom &
  !       ( LogE,    LogT,    LogEta, &
  !         LogEs_T, LogTs_T, LogEtas_T, &
  !         OS_NES(1,iH1), NES_T(:,:,:,:,iH1), opH1_CPU_1 )
  !tEnd = get_wtime()
  tCPU = 0.0d0

  tBegin = get_wtime()
  CALL LogInterpolateSingleVariable_2D2D_Custom &
         ( LogE,    LogT,    LogEta, &
           LogEs_T, LogTs_T, LogEtas_T, &
           OS_NES(1,iH1), NES_T(:,:,:,:,iH1), opH1_GPU_1 )
  tEnd = get_wtime()
  tGPU = tEnd - tBegin

  WRITE(*,*)
  WRITE(*,'(A4,A64,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_2D2D_Custom (CPU): ', tCPU
  WRITE(*,'(A4,A64,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_2D2D_Custom (GPU): ', tGPU
  WRITE(*,*)

  tBegin = get_wtime()
  DO iP = 1, nPointsX

    CALL LogInterpolateSingleVariable_2D2D_Custom_Point &
           ( LogE,    LogT(iP), LogEta(iP), &
             LogEs_T, LogTs_T,  LogEtas_T,   &
             OS_NES(1,iH1), NES_T(:,:,:,:,iH1), opH1_CPU_2(:,:,iP) )

  END DO
  tEnd = get_wtime()
  tCPU = tEnd - tBegin
  opH1_CPU_1 = opH1_CPU_2

  tBegin = get_wtime()
#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined (WEAKLIB_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR &
  !$ACC PRESENT( LogE, LogT, LogEta, LogEs_T, LogTs_T, LogEtas_T, &
  !$ACC          OS_NES, NES_T, opH1_GPU_2 )
#endif
  DO iP = 1, nPointsX

    CALL LogInterpolateSingleVariable_2D2D_Custom_Point &
           ( LogE,    LogT(iP), LogEta(iP), &
             LogEs_T, LogTs_T,  LogEtas_T,   &
             OS_NES(1,iH1), NES_T(:,:,:,:,iH1), opH1_GPU_2(:,:,iP) )

  END DO
  tEnd = get_wtime()
  tGPU = tEnd - tBegin

  WRITE(*,*)
  WRITE(*,'(A4,A64,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_2D2D_Custom_Point (CPU): ', tCPU
  WRITE(*,'(A4,A64,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_2D2D_Custom_Point (GPU): ', tGPU
  WRITE(*,*)

  !tBegin = get_wtime()
  !CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
  !       ( LogT,    LogEta, &
  !         LogTs_T, LogEtas_T, &
  !         OS_NES(1,iH1), NES_AT(:,:,:,:,iH1), opH1_CPU_3 )
  !tEnd = get_wtime()
  tCPU = 0.0d0

  tBegin = get_wtime()
  CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned &
         ( LogT,    LogEta, &
           LogTs_T, LogEtas_T, &
           OS_NES(1,iH1), NES_AT(:,:,:,:,iH1), opH1_GPU_3 )
  tEnd = get_wtime()
  tGPU = tEnd - tBegin

  WRITE(*,*)
  WRITE(*,'(A4,A64,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_2D2D_Custom_Aligned (CPU): ', tCPU
  WRITE(*,'(A4,A64,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_2D2D_Custom_Aligned (GPU): ', tGPU
  WRITE(*,*)

  tBegin = get_wtime()
  DO iP = 1, nPointsX

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned_Point &
           ( LogT(iP), LogEta(iP), &
             LogTs_T,  LogEtas_T, &
             OS_NES(1,iH1), NES_AT(:,:,:,:,iH1), opH1_CPU_4(:,:,iP) )

  END DO
  tEnd = get_wtime()
  tCPU = tEnd - tBegin
  opH1_CPU_3 = opH1_CPU_4

  tBegin = get_wtime()
#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD
#elif defined (WEAKLIB_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR &
  !$ACC PRESENT( LogT, LogEta, LogTs_T, LogEtas_T, &
  !$ACC          OS_NES, NES_AT, opH1_GPU_4 )
#endif
  DO iP = 1, nPointsX

    CALL LogInterpolateSingleVariable_2D2D_Custom_Aligned_Point &
           ( LogT(iP), LogEta(iP), &
             LogTs_T,  LogEtas_T, &
             OS_NES(1,iH1), NES_AT(:,:,:,:,iH1), opH1_GPU_4(:,:,iP) )

  END DO
  tEnd = get_wtime()
  tGPU = tEnd - tBegin

  WRITE(*,*)
  WRITE(*,'(A4,A64,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_2D2D_Custom_Aligned_Point (CPU): ', tCPU
  WRITE(*,'(A4,A64,ES10.4E2)') &
    '', 'LogInterpolateSingleVariable_2D2D_Custom_Aligned_Point (GPU): ', tGPU
  WRITE(*,*)

  WRITE(*,*) 'NES (CPU)'
  WRITE(*,*) MINVAL( opH1_CPU_1 ), MAXVAL( opH1_CPU_1 )
  WRITE(*,*) MINVAL( opH1_CPU_2 ), MAXVAL( opH1_CPU_2 )
  WRITE(*,*) MINVAL( opH1_CPU_3 ), MAXVAL( opH1_CPU_3 )
  WRITE(*,*) MINVAL( opH1_CPU_4 ), MAXVAL( opH1_CPU_4 )
  WRITE(*,*) MAXVAL( ABS( opH1_CPU_1 - opH1_CPU_2 ) )
  WRITE(*,*) MAXVAL( ABS( opH1_CPU_1 - opH1_CPU_3 ) )
  WRITE(*,*) MAXVAL( ABS( opH1_CPU_1 - opH1_CPU_4 ) )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE FROM( opH1_GPU_1, opH1_GPU_2, opH1_GPU_3, opH1_GPU_4 )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE HOST( opH1_GPU_1, opH1_GPU_2, opH1_GPU_3, opH1_GPU_4 )
#endif

  WRITE(*,*) 'NES (GPU)'
  WRITE(*,*) MINVAL( opH1_GPU_1 ), MAXVAL( opH1_GPU_1 )
  WRITE(*,*) MINVAL( opH1_GPU_2 ), MAXVAL( opH1_GPU_2 )
  WRITE(*,*) MINVAL( opH1_GPU_3 ), MAXVAL( opH1_GPU_3 )
  WRITE(*,*) MINVAL( opH1_GPU_4 ), MAXVAL( opH1_GPU_4 )
  WRITE(*,*) MAXVAL( ABS( opH1_GPU_1 - opH1_GPU_2 ) )
  WRITE(*,*) MAXVAL( ABS( opH1_GPU_1 - opH1_GPU_3 ) )
  WRITE(*,*) MAXVAL( ABS( opH1_GPU_1 - opH1_GPU_4 ) )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET EXIT DATA &
  !$OMP MAP( release: LogEs_T, LogDs_T, LogTs_T, Ys_T, LogEtas_T, &
  !$OMP               LogE, LogD, LogT, Y, LogEta, &
  !$OMP               OS_EmAb, EmAb_T, OS_NES, NES_T, NES_AT, &
  !$OMP               opEC_GPU_1, opEC_GPU_2, opEC_GPU_3, &
  !$OMP               opH1_GPU_1, opH1_GPU_2, opH1_GPU_3, opH1_GPU_4 )
#elif defined (WEAKLIB_OACC)
  !$ACC EXIT DATA &
  !$ACC DELETE( LogEs_T, LogDs_T, LogTs_T, Ys_T, LogEtas_T, &
  !$ACC         LogE, LogD, LogT, Y, LogEta, &
  !$ACC         OS_EmAb, EmAb_T, OS_NES, NES_T, NES_AT, &
  !$ACC         opEC_GPU_1, opEC_GPU_2, opEC_GPU_3, &
  !$ACC         opH1_GPU_1, opH1_GPU_2, opH1_GPU_3, opH1_GPU_4 )
#endif

  END ASSOCIATE

#if defined(WEAKLIB_OMP_OL)
#elif defined(WEAKLIB_OACC)
  !!$ACC SHUTDOWN
#endif

  CALL MPI_FINALIZE( ierr )

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
