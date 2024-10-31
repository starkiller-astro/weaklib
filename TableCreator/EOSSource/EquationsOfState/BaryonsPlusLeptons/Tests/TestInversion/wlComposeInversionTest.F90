PROGRAM wlComposeInversionTest

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
  USE wlLeptonEOSModule, ONLY: &
    HelmholtzEOSType, MuonEOSType
  USE wlElectronEOS, ONLY: &
    FullHelmEOS, MinimalHelmEOS_rt, ElectronStateType
  USE wlMuonEOS, ONLY: &
    FullMuonEOS, MuonStateType
  USE wlHelmMuonIOModuleHDF, ONLY: &
    ReadHelmholtzTableHDF, ReadMuonTableHDF
  USE wlExtPhysicalConstantsModule, ONLY: &
    kmev, rmu, kmev_inv, ergmev, me, mmu, cvel
    
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
    EOSBaryonTable
  TYPE(HelmholtzEOSType) :: HelmholtzTable
  TYPE(ElectronStateType) :: ElectronState
  TYPE(MuonEOSType) :: MuonTable
  TYPE(MuonStateType) :: MuonState
  REAL(DP), DIMENSION(:), ALLOCATABLE :: &
    Ds_T, Ts_T, Ys_T
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: &
    Ps_T, Es_T, Ss_T
  REAL(DP) :: &
    OS_P, OS_E, OS_S
  INTEGER :: &
    iD_T, iT_T, iY_T, iP_T, iS_T, iE_T
  
  INTEGER :: iRho, iTemp, iYp
  REAL(DP) :: Ps_T_Helm, Es_T_Helm, Ss_T_Helm
  REAL(DP) :: Ps_T_muon, Es_T_muon, Ss_T_muon, Ymu
  LOGICAL :: OnlyOneTable
  
  CHARACTER(len=128) :: BaryonPlusHelmTableName, FullTableName
  
  CALL MPI_INIT( ierr )

! #if defined(WEAKLIB_OMP_OL)
! #elif defined(WEAKLIB_OACC)
  ! !!$ACC INIT
! #endif

  OnlyOneTable = .true.
  
  BaryonPlusHelmTableName = '../../../../../../../../weaklib_tables/BaryonsPlusHelmPlusMuonsEOS_interpolated.h5'
  BaryonPlusHelmTableName = '../../../../../../../../weaklib_tables/BaryonsPlusHelmPlusMuonsEOS.h5'
  BaryonPlusHelmTableName = "/mnt/c/Users/Lucab/GDrive/Research/Muons_project/weaklib_tables/FullEOS.h5"

  CALL InitializeHDF( )
  
  IF (OnlyOneTable) THEN
    CALL ReadEquationOfStateTableHDF( EOSBaryonTable, BaryonPlusHelmTableName )
  ELSE
    ! read in helmholtz table
    CALL ReadHelmholtzTableHDF( HelmholtzTable, BaryonPlusHelmTableName )

    ! read in helmholtz table
    CALL ReadMuonTableHDF( MuonTable, BaryonPlusHelmTableName )

    ! read in baryon table -------------------------------
    CALL ReadEquationOfStateTableHDF( EOSBaryonTable, BaryonPlusHelmTableName )
  ENDIF
  
  CALL FinalizeHDF( )

  iD_T = EOSBaryonTable % TS % Indices % iRho
  iT_T = EOSBaryonTable % TS % Indices % iT
  iY_T = EOSBaryonTable % TS % Indices % iYe

  ALLOCATE( Ds_T(EOSBaryonTable % TS % nPoints(iD_T)) )
  Ds_T = EOSBaryonTable % TS % States(iD_T) % Values

  ALLOCATE( Ts_T(EOSBaryonTable % TS % nPoints(iT_T)) )
  Ts_T = EOSBaryonTable % TS % States(iT_T) % Values

  ALLOCATE( Ys_T(EOSBaryonTable % TS % nPoints(iY_T)) )
  Ys_T = EOSBaryonTable % TS % States(iY_T) % Values

  iP_T = EOSBaryonTable % DV % Indices % iPressure
  iS_T = EOSBaryonTable % DV % Indices % iEntropyPerBaryon
  iE_T = EOSBaryonTable % DV % Indices % iInternalEnergyDensity

  OS_P = EOSBaryonTable % DV % Offsets(iP_T)
  OS_S = EOSBaryonTable % DV % Offsets(iS_T)
  OS_E = EOSBaryonTable % DV % Offsets(iE_T)

  ALLOCATE &
    ( Ps_T(1:EOSBaryonTable % DV % nPoints(1), &
           1:EOSBaryonTable % DV % nPoints(2), &
           1:EOSBaryonTable % DV % nPoints(3)) )
  ALLOCATE &
    ( Ss_T(1:EOSBaryonTable % DV % nPoints(1), &
           1:EOSBaryonTable % DV % nPoints(2), &
           1:EOSBaryonTable % DV % nPoints(3)) )
  ALLOCATE &
    ( Es_T(1:EOSBaryonTable % DV % nPoints(1), &
           1:EOSBaryonTable % DV % nPoints(2), &
           1:EOSBaryonTable % DV % nPoints(3)) )

  Ps_T = EOSBaryonTable % DV % Variables(iP_T ) % Values
  Ss_T = EOSBaryonTable % DV % Variables(iS_T ) % Values
  Es_T = EOSBaryonTable % DV % Variables(iE_T ) % Values

  IF (.not. OnlyOneTable) THEN
    DO iRho=1,EOSBaryonTable % DV % nPoints(1)
      DO iTemp=1,EOSBaryonTable % DV % nPoints(2)
        DO iYp=1,EOSBaryonTable % DV % nPoints(3)
  
          Ymu = Ys_T(iYp) / 5.0d0

          ! Now add electron component
          ! Initialize temperature, density, yp, Zbar and Abar
          ElectronState % t = Ts_T(iTemp)
          ElectronState % rho = Ds_T(iRho)
          ElectronState % abar = 1.0d0 ! these are only used for ion contribution
          ElectronState % zbar = 1.0d0 ! these are only used for ion contribution
          ElectronState % y_e = Ys_T(iYp) - Ymu
          
          ! calculate electron quantities
          CALL FullHelmEOS(1, HelmholtzTable, ElectronState, .false., .false.)
          CALL MinimalHelmEOS_rt(HelmholtzTable, ElectronState)

          Es_T_helm = ElectronState % e + me / rmu * ergmev * ElectronState % y_e ! add back mass to internal energy!
          Ps_T_helm = ElectronState % p
          Ss_T_helm = ElectronState % s

          ! calculate muon quantities 
          MuonState % t = Ts_T(iTemp)
          MuonState % rhoymu = Ds_T(iRho) * Ymu
          
          CALL FullMuonEOS(MuonTable, MuonState)

          Es_T_muon = MuonState % e + mmu / rmu * ergmev * Ymu ! add back mass to internal energy!
          Ps_T_muon = MuonState % p
          Ss_T_muon = MuonState % s
                  
          Es_T(iRho,iTemp,iYp) = LOG10( 10.0d0**( Es_T(iRho,iTemp,iYp) ) + Es_T_helm + Es_T_muon)
          Ps_T(iRho,iTemp,iYp) = LOG10( 10.0d0**( Ps_T(iRho,iTemp,iYp) ) + Ps_T_helm + Ps_T_muon)
          Ss_T(iRho,iTemp,iYp) = LOG10( 10.0d0**( Ss_T(iRho,iTemp,iYp) ) + Ss_T_helm + Ss_T_muon)

          IF (Es_T(iRho,iTemp,iYp) .NE. Es_T(iRho,iTemp,iYp)) THEN
            WRITE(*,*) 'E', iRho,iTemp,iYp, Es_T(iRho,iTemp,iYp), Es_T_helm, Es_T_muon
            STOP
          ENDIF
          IF (Ps_T(iRho,iTemp,iYp) .NE. Ps_T(iRho,iTemp,iYp)) THEN
            WRITE(*,*) 'P', iRho,iTemp,iYp, Ps_T(iRho,iTemp,iYp), Ps_T_helm, Ps_T_muon
            STOP
          ENDIF
          IF (Ss_T(iRho,iTemp,iYp) .NE. Ss_T(iRho,iTemp,iYp)) THEN
            WRITE(*,*) 'S', iRho,iTemp,iYp, Ss_T(iRho,iTemp,iYp), Ss_T_helm, Ss_T_muon
            STOP
          ENDIF
          
          ENDDO
        ENDDO
      ENDDO
      
  ENDIF

  CALL InitializeEOSInversion &
         ( Ds_T, Ts_T, Ys_T, &
           10.0d0**( Es_T ) - OS_E, &
           10.0d0**( Ps_T ) - OS_P, &
           10.0d0**( Ss_T ) - OS_S, &
           Verbose_Option = .TRUE. )

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET ENTER DATA &
  ! !$OMP MAP( to: Ds_T, Ts_T, Ys_T, Es_T, Ps_T, Ss_T, OS_E, OS_P, OS_S ) &
  ! !$OMP MAP( alloc: D, T, Y, P, E, S, T_P, T_E, T_S, Error_P, Error_E, Error_S )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC ENTER DATA &
  ! !$ACC COPYIN( Ds_T, Ts_T, Ys_T, Es_T, Ps_T, Ss_T, OS_E, OS_P, OS_S ) &
  ! !$ACC CREATE( D, T, Y, P, E, S, T_P, T_E, T_S, Error_P, Error_E, Error_S )
! #endif

  WRITE(*,*)
  WRITE(*,'(A4,A10,I10.10)') '', 'nPoints = ', nPoints

  ! --- Random Sampling of EOSBaryonTable Table ---

  n_rndm = nPoints

  ! --- Initialize Density Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_D )

  associate &
    ( minD => EOSBaryonTable % TS % MinValues(iD), &
      maxD => EOSBaryonTable % TS % MaxValues(iD) )

  D(:) = 10**( LOG10(minD) + ( LOG10(maxD) - LOG10(minD) ) * rndm_D )

  PRINT*, "Min/Max D = ", MINVAL( D ), MAXVAL( D )
  PRINT*, "            ", minD, maxD

  end associate

  ! --- Initialize Temperature Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_T )

  associate &
    ( minT => EOSBaryonTable % TS % MinValues(iT), &
      maxT => EOSBaryonTable % TS % MaxValues(iT) )

  T(:) = 10**( LOG10(minT) + ( LOG10(maxT) - LOG10(minT) ) * rndm_T )

  PRINT*, "Min/Max T = ", MINVAL( T ), MAXVAL( T )
  PRINT*, "            ", minT, maxT

  end associate

  ! T(:) = 10**( LOG10(3.0d0*kmev_inv) + ( LOG10(30.0d0*kmev_inv) - LOG10(3.0d0*kmev_inv) ) * rndm_T )
  PRINT*, "Min/Max T = ", MINVAL( T ), MAXVAL( T )

  ! --- Initialize Electron Fraction Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_Y )

  associate &
    ( minY => EOSBaryonTable % TS % MinValues(iY), &
      maxY => EOSBaryonTable % TS % MaxValues(iY) )

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
  WRITE(*,*) "Min/Max E = ", MINVAL( Es_T ), MAXVAL( Es_T )
  WRITE(*,*) "Min/Max P = ", MINVAL( P ), MAXVAL( P )
  WRITE(*,*) "Min/Max P = ", MINVAL( Ps_T ), MAXVAL( Ps_T )
  WRITE(*,*) "Min/Max S = ", MINVAL( S ), MAXVAL( S )
  WRITE(*,*) "Min/Max S = ", MINVAL( Ss_T ), MAXVAL( Ss_T )

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET UPDATE TO( D, T, Y, P, E, S )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC UPDATE DEVICE( D, T, Y, P, E, S )
! #endif

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
    ( minT => 1.0001 * EOSBaryonTable % TS % MinValues(iT), &
      maxT => 0.9999 * EOSBaryonTable % TS % MaxValues(iT) )

  T_E = MIN( MAX( T + dT, minT ), maxT )

  end associate

  Error_E = 0

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET UPDATE TO( T_E, Error_E )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC UPDATE DEVICE( T_E, Error_E )
! #endif

  CALL CPU_TIME( tBegin )

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  ! !$OMP PRIVATE( T_Guess )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC PARALLEL LOOP GANG VECTOR &
  ! !$ACC PRIVATE( T_Guess ) &
  ! !$ACC PRESENT( D, E, Y, Ds_T, Ts_T, Ys_T, Es_T, OS_E, T_E, Error_E )
! #elif defined (WEAKLIB_OMP)
  ! !$OMP PARALLEL DO &
  ! !$OMP PRIVATE( T_Guess )
! #endif
  DO iP = 1, nPoints
    T_Guess = T_E(iP)
    CALL ComputeTemperatureWith_DEY_Single_Guess &
           ( D(iP), E(iP), Y(iP), Ds_T, Ts_T, Ys_T, Es_T, OS_E, T_E(iP), T_Guess, &
             Error_E(iP) )
  END DO

  CALL CPU_TIME( tEnd )
  tGPU = tEnd - tBegin

  PRINT*
  PRINT*, "ComputeTemperatureWith_DEYp (Good Guess):"
  PRINT*, "GPU_TIME = ", tGPU

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
  STOP
  
  ! -------------------------------------------------------------------

  associate &
    ( minT => 1.0001 * EOSBaryonTable % TS % MinValues(iT), &
      maxT => 0.9999 * EOSBaryonTable % TS % MaxValues(iT) )

  Amp = 0.1_dp * ( maxT - minT ) ! --- Perturbation Amplitude

  PRINT*
  PRINT*, "Perturbation Amplitude: ", Amp

  dT = Amp * 2.0_dp * ( rndm_T - 0.5_dp ) * T

  T_E = MIN( MAX( T + dT, minT ), maxT )

  end associate

  Error_E = 0

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET UPDATE TO( T_E, Error_E )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC UPDATE DEVICE( T_E, Error_E )
! #endif

  CALL CPU_TIME( tBegin )

  CALL ComputeTemperatureWith_DEY_Many &
         ( D, E, Y, Ds_T, Ts_T, Ys_T, Es_T, OS_E, T_E, &
           UseInitialGuess_Option = .TRUE., &
           Error_Option = Error_E )

  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  CALL CPU_TIME( tBegin )

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  ! !$OMP PRIVATE( T_Guess )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC PARALLEL LOOP GANG VECTOR &
  ! !$ACC PRIVATE( T_Guess ) &
  ! !$ACC PRESENT( D, E, Y, Ds_T, Ts_T, Ys_T, Es_T, OS_E, T_E, Error_E )
! #elif defined (WEAKLIB_OMP)
  ! !$OMP PARALLEL DO &
  ! !$OMP PRIVATE( T_Guess )
! #endif
  DO iP = 1, nPoints
    T_Guess = T_E(iP)
    CALL ComputeTemperatureWith_DEY_Single_Guess &
           ( D(iP), E(iP), Y(iP), Ds_T, Ts_T, Ys_T, Es_T, OS_E, T_E(iP), T_Guess, &
             Error_E(iP) )
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

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET UPDATE FROM( T_E, Error_E )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC UPDATE HOST( T_E, Error_E )
! #endif
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

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET UPDATE TO( T_E, Error_E )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC UPDATE DEVICE( T_E, Error_E )
! #endif

  CALL CPU_TIME( tBegin )

  CALL ComputeTemperatureWith_DEY_Many &
         ( D, E, Y, Ds_T, Ts_T, Ys_T, Es_T, OS_E, T_E, &
           UseInitialGuess_Option = .FALSE., &
           Error_Option = Error_E )

  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  CALL CPU_TIME( tBegin )

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  ! !$OMP PRIVATE( T_Guess )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC PARALLEL LOOP GANG VECTOR &
  ! !$ACC PRIVATE( T_Guess ) &
  ! !$ACC PRESENT( D, E, Y, Ds_T, Ts_T, Ys_T, Es_T, OS_E, T_E, Error_E )
! #elif defined (WEAKLIB_OMP)
  ! !$OMP PARALLEL DO &
  ! !$OMP PRIVATE( T_Guess )
! #endif
  DO iP = 1, nPoints
    T_Guess = T_E(iP)
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
  Error = ABS( T - T_E ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_E(iMaxError), Y(iMaxError), Ds_T, Ts_T, Ys_T, OS_E, Es_T, E_T )
  PRINT*, "CPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "CPU  Error(E) = ", ABS( E(iMaxError) - E_T ) / E(iMaxError)
  PRINT*, "CPU MINVAL(T) = ", MINVAL( T_E )

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET UPDATE FROM( T_E, Error_E )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC UPDATE HOST( T_E, Error_E )
! #endif
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
    ( minT => 1.0001 * EOSBaryonTable % TS % MinValues(iT), &
      maxT => 0.9999 * EOSBaryonTable % TS % MaxValues(iT) )

  T_S = MIN( MAX( T + dT, minT ), maxT )

  end associate

  Error_S = 0

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET UPDATE TO( T_S, Error_S )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC UPDATE DEVICE( T_S, Error_S )
! #endif

  CALL CPU_TIME( tBegin )

  CALL ComputeTemperatureWith_DSY_Many &
         ( D, S, Y, Ds_T, Ts_T, Ys_T, Ss_T, OS_S, T_S, &
           UseInitialGuess_Option = .TRUE., &
           Error_Option = Error_S )

  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  CALL CPU_TIME( tBegin )

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  ! !$OMP PRIVATE( T_Guess )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC PARALLEL LOOP GANG VECTOR &
  ! !$ACC PRIVATE( T_Guess ) &
  ! !$ACC PRESENT( D, S, Y, Ds_T, Ts_T, Ys_T, Ss_T, OS_S, T_S, Error_S )
! #elif defined (WEAKLIB_OMP)
  ! !$OMP PARALLEL DO &
  ! !$OMP PRIVATE( T_Guess )
! #endif
  DO iP = 1, nPoints
    T_Guess = T_S(iP)
    CALL ComputeTemperatureWith_DSY_Single_Guess &
           ( D(iP), S(iP), Y(iP), Ds_T, Ts_T, Ys_T, Ss_T, OS_S, T_S(iP), T_Guess, &
             Error_S(iP) )
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

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET UPDATE FROM( T_S, Error_S )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC UPDATE HOST( T_S, Error_S )
! #endif
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
    ( minT => 1.0001 * EOSBaryonTable % TS % MinValues(iT), &
      maxT => 0.9999 * EOSBaryonTable % TS % MaxValues(iT) )

  Amp = 0.1_dp * ( maxT - minT )
  dT = Amp * 2.0_dp * ( rndm_T - 0.5_dp ) * T

  T_S = MIN( MAX( T + dT, minT ), maxT )

  end associate

  Error_S = 0

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET UPDATE TO( T_S, Error_S )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC UPDATE DEVICE( T_S, Error_S )
! #endif

  CALL CPU_TIME( tBegin )

  CALL ComputeTemperatureWith_DSY_Many &
         ( D, S, Y, Ds_T, Ts_T, Ys_T, Ss_T, OS_S, T_S, &
           UseInitialGuess_Option = .TRUE., &
           Error_Option = Error_S )

  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  CALL CPU_TIME( tBegin )

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  ! !$OMP PRIVATE( T_Guess )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC PARALLEL LOOP GANG VECTOR &
  ! !$ACC PRIVATE( T_Guess ) &
  ! !$ACC PRESENT( D, S, Y, Ds_T, Ts_T, Ys_T, Ss_T, OS_S, T_S, Error_S )
! #elif defined (WEAKLIB_OMP)
  ! !$OMP PARALLEL DO &
  ! !$OMP PRIVATE( T_Guess )
! #endif
  DO iP = 1, nPoints
    T_Guess = T_S(iP)
    CALL ComputeTemperatureWith_DSY_Single_Guess &
           ( D(iP), S(iP), Y(iP), Ds_T, Ts_T, Ys_T, Ss_T, OS_S, T_S(iP), T_Guess, &
             Error_S(iP) )
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

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET UPDATE FROM( T_S, Error_S )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC UPDATE HOST( T_S, Error_S )
! #endif
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

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET UPDATE TO( T_S, Error_S )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC UPDATE DEVICE( T_S, Error_S )
! #endif

  CALL CPU_TIME( tBegin )

  CALL ComputeTemperatureWith_DSY_Many &
         ( D, S, Y, Ds_T, Ts_T, Ys_T, Ss_T, OS_S, T_S, &
           UseInitialGuess_Option = .FALSE., &
           Error_Option = Error_S )

  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  CALL CPU_TIME( tBegin )

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  ! !$OMP PRIVATE( T_Guess )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC PARALLEL LOOP GANG VECTOR &
  ! !$ACC PRIVATE( T_Guess ) &
  ! !$ACC PRESENT( D, S, Y, Ds_T, Ts_T, Ys_T, Ss_T, OS_S, T_S, Error_S )
! #elif defined (WEAKLIB_OMP)
  ! !$OMP PARALLEL DO &
  ! !$OMP PRIVATE( T_Guess )
! #endif
  DO iP = 1, nPoints
    T_Guess = T_S(iP)
    CALL ComputeTemperatureWith_DSY_Single_NoGuess &
           ( D(iP), S(iP), Y(iP), Ds_T, Ts_T, Ys_T, Ss_T, OS_S, T_S(iP), &
             Error_S(iP) )
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

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET UPDATE FROM( T_S, Error_S )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC UPDATE HOST( T_S, Error_S )
! #endif
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
    ( minT => 1.0001 * EOSBaryonTable % TS % MinValues(iT), &
      maxT => 0.9999 * EOSBaryonTable % TS % MaxValues(iT) )

  T_P = MIN( MAX( T + dT, minT ), maxT )

  end associate

  Error_P = 0

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET UPDATE TO( T_P, Error_P )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC UPDATE DEVICE( T_P, Error_P )
! #endif

  CALL CPU_TIME( tBegin )

  CALL ComputeTemperatureWith_DPY_Many &
         ( D, P, Y, Ds_T, Ts_T, Ys_T, Ps_T, OS_P, T_P, &
           UseInitialGuess_Option = .TRUE., &
           Error_Option = Error_P )

  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  CALL CPU_TIME( tBegin )

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  ! !$OMP PRIVATE( T_Guess )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC PARALLEL LOOP GANG VECTOR &
  ! !$ACC PRIVATE( T_Guess ) &
  ! !$ACC PRESENT( D, P, Y, Ds_T, Ts_T, Ys_T, Ps_T, OS_P, T_P, Error_P )
! #elif defined(WEAKLIB_OMP_OL)
  ! !$OMP PARALLEL DO &
  ! !$OMP PRIVATE( T_Guess )
! #endif
  DO iP = 1, nPoints
    T_Guess = T_P(iP)
    CALL ComputeTemperatureWith_DPY_Single_Guess &
           ( D(iP), P(iP), Y(iP), Ds_T, Ts_T, Ys_T, Ps_T, OS_P, T_P(iP), T_Guess, &
             Error_P(iP) )
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

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET UPDATE FROM( T_P, Error_P )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC UPDATE HOST( T_P, Error_P )
! #endif
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
    ( minT => 1.0001 * EOSBaryonTable % TS % MinValues(iT), &
      maxT => 0.9999 * EOSBaryonTable % TS % MaxValues(iT) )

  Amp = 0.1_dp * ( maxT - minT ) ! --- Perturbation Amplitude
  dT = Amp * 2.0_dp * ( rndm_T - 0.5_dp ) * T

  T_P = MIN( MAX( T + dT, minT ), maxT )

  end associate

  Error_P = 0

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET UPDATE TO( T_P, Error_P )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC UPDATE DEVICE( T_P, Error_P )
! #endif

  CALL CPU_TIME( tBegin )

  CALL ComputeTemperatureWith_DPY_Many &
         ( D, P, Y, Ds_T, Ts_T, Ys_T, Ps_T, OS_P, T_P, &
           UseInitialGuess_Option = .TRUE., &
           Error_Option = Error_P )

  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  CALL CPU_TIME( tBegin )

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  ! !$OMP PRIVATE( T_Guess )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC PARALLEL LOOP GANG VECTOR &
  ! !$ACC PRIVATE( T_Guess ) &
  ! !$ACC PRESENT( D, P, Y, Ds_T, Ts_T, Ys_T, Ps_T, OS_P, T_P, Error_P )
! #elif defined(WEAKLIB_OMP_OL)
  ! !$OMP PARALLEL DO &
  ! !$OMP PRIVATE( T_Guess )
! #endif
  DO iP = 1, nPoints
    T_Guess = T_P(iP)
    CALL ComputeTemperatureWith_DPY_Single_Guess &
           ( D(iP), P(iP), Y(iP), Ds_T, Ts_T, Ys_T, Ps_T, OS_P, T_P(iP), T_Guess, &
             Error_P(iP) )
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

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET UPDATE FROM( T_P, Error_P )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC UPDATE HOST( T_P, Error_P )
! #endif
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

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET UPDATE TO( T_P, Error_P )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC UPDATE DEVICE( T_P, Error_P )
! #endif

  CALL CPU_TIME( tBegin )

  CALL ComputeTemperatureWith_DPY_Many &
         ( D, P, Y, Ds_T, Ts_T, Ys_T, Ps_T, OS_P, T_P, &
           UseInitialGuess_Option = .FALSE., &
           Error_Option = Error_P )

  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  CALL CPU_TIME( tBegin )

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  ! !$OMP PRIVATE( T_Guess )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC PARALLEL LOOP GANG VECTOR &
  ! !$ACC PRIVATE( T_Guess ) &
  ! !$ACC PRESENT( D, P, Y, Ds_T, Ts_T, Ys_T, Ps_T, OS_P, T_P, Error_P )
! #elif defined (WEAKLIB_OMP)
  ! !$OMP PARALLEL DO &
  ! !$OMP PRIVATE( T_Guess )
! #endif
  DO iP = 1, nPoints
    T_Guess = T_P(iP)
    CALL ComputeTemperatureWith_DPY_Single_NoGuess &
           ( D(iP), P(iP), Y(iP), Ds_T, Ts_T, Ys_T, Ps_T, OS_P, T_P(iP), &
             Error_P(iP) )
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

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET UPDATE FROM( T_P, Error_P )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC UPDATE HOST( T_P, Error_P )
! #endif
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

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET EXIT DATA &
  ! !$OMP MAP( release: D, E, P, S, Y, Ds_T, Ts_T, Ys_T, Es_T, Ps_T, Ss_T, OS_E, OS_P, OS_S, T_E, Error_E )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC EXIT DATA &
  ! !$ACC DELETE( D, E, P, S, Y, Ds_T, Ts_T, Ys_T, Es_T, Ps_T, Ss_T, OS_E, OS_P, OS_S, T_E, Error_E )
! #endif

  DEALLOCATE( Ds_T, Ts_T, Ys_T, Ps_T, Ss_T, Es_T )

! #if defined(WEAKLIB_OMP_OL)
! #elif defined(WEAKLIB_OACC)
  ! !!$ACC SHUTDOWN
! #endif

  CALL MPI_FINALIZE( ierr )

END PROGRAM wlComposeInversionTest
