PROGRAM wlComposeNewInversionTest

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
  USE wlInterpolationUtilitiesModule, ONLY: &
    Index1D_Lin, &
    Index1D_Log
#ifdef INTERPOLATION_SPLIT_TABLE_COMBINED
  USE wlEOSComponentsCombinedInversionModule
#else
  USE wlEOSComponentsSeparateInversionModule
#endif
  USE wlLeptonEOSModule, ONLY: &
    HelmholtzTableType, MuonEOSType
  USE wlElectronPhotonEOS, ONLY: &
    ElectronPhotonEOS, ElectronPhotonStateType
  USE wlMuonEOS, ONLY: &
    FullMuonEOS, MuonStateType
  USE wlHelmMuonIOModuleHDF, ONLY: &
    ReadHelmholtzTableHDF, ReadMuonTableHDF
  USE wlExtPhysicalConstantsModule, ONLY: &
    kmev, rmu, kmev_inv, ergmev, cvel
    
  IMPLICIT NONE

  INCLUDE 'mpif.h'

  INTEGER :: &
    n_rndm, &
    iP, &
    iMaxError, &
    ierr
  INTEGER, PARAMETER :: &
    nPoints = 2**16
  INTEGER, DIMENSION(nPoints) :: &
    Error_E, &
    Error_P, &
    Error_S
  REAL(dp) :: &
    tBegin, &
    tEnd, &
    tCPU, &
    T_Guess, &
    Amp, &
    E_T, S_T, P_T
  REAL(dp), DIMENSION(nPoints) :: &
    D, T, Yp, Ye, Ym, E, P, S, T_E, T_P, T_S, dT, &
    E_int_error, P_int_error, S_int_error, &
    rndm_D, &
    rndm_T, &
    rndm_Yp, &
    Error
  TYPE(EquationOfStateTableType) :: &
    EOSBaryonTable
  TYPE(HelmholtzTableType) :: HelmholtzTable
  TYPE(ElectronPhotonStateType) :: ElectronPhotonState
  TYPE(MuonEOSType) :: MuonTable
  TYPE(MuonStateType) :: MuonState
  REAL(DP), DIMENSION(:), ALLOCATABLE :: &
    Ds_bary, Ts_bary, Yps_bary, &
    Ds_full, Ts_full, Yps_full
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: &
    Ps_bary, Es_bary, Ss_bary, &
    Ps_full, Es_full, Ss_full
  REAL(DP) :: &
    OS_P, OS_E, OS_S
  INTEGER :: SizeDs, SizeTs, SizeYps, &
    iD_bary, iT_bary, iYp_bary, iP_bary, iS_bary, iE_bary
  
  INTEGER  :: iD, iT, iYp, iwrong
  REAL(DP) :: LocalOffset
  REAL(DP) :: P_bary, E_bary, S_bary
  REAL(DP) :: Ps_Helm, Es_Helm, Ss_Helm
  REAL(DP) :: Ymu_temp, Ps_muon, Es_muon, Ss_muon, Yp_over_Ymu

  CHARACTER(len=128) :: BaryonPlusHelmTableName
    
  CALL MPI_INIT( ierr )

#if defined(WEAKLIB_OMP_OL)
#elif defined(WEAKLIB_OACC)
  !!$ACC INIT
#endif
  
  BaryonPlusHelmTableName = 'BaryonsPlusHelmPlusMuonsEOS_interpolated.h5'
  !BaryonPlusHelmTableName = 'BaryonsPlusHelmPlusMuonsEOS.h5'
  Yp_over_Ymu = 0.0d0
  
  CALL InitializeHDF( )
  
  ! read in helmholtz table
  CALL ReadHelmholtzTableHDF( HelmholtzTable, BaryonPlusHelmTableName )

  ! read in muon table
  CALL ReadMuonTableHDF( MuonTable, BaryonPlusHelmTableName )

  ! read in baryon table -------------------------------
  CALL ReadEquationOfStateTableHDF( EOSBaryonTable, BaryonPlusHelmTableName )
  
  CALL FinalizeHDF( )

  iD_bary  = EOSBaryonTable % TS % Indices % iRho
  iT_bary  = EOSBaryonTable % TS % Indices % iT
  iYp_bary = EOSBaryonTable % TS % Indices % iYe

  SizeDs  = EOSBaryonTable % TS % nPoints(iD_bary )
  SizeTs  = EOSBaryonTable % TS % nPoints(iT_bary )
  SizeYps = EOSBaryonTable % TS % nPoints(iYp_bary)

  ALLOCATE( Ds_bary(SizeDs) )
  Ds_bary = EOSBaryonTable % TS % States(iD_bary) % Values

  ALLOCATE( Ts_bary(SizeTs) )
  Ts_bary = EOSBaryonTable % TS % States(iT_bary) % Values

  ALLOCATE( Yps_bary(SizeYps) )
  Yps_bary = EOSBaryonTable % TS % States(iYp_bary) % Values

  iP_bary = EOSBaryonTable % DV % Indices % iPressure
  iS_bary = EOSBaryonTable % DV % Indices % iEntropyPerBaryon
  iE_bary = EOSBaryonTable % DV % Indices % iInternalEnergyDensitY

  OS_P = EOSBaryonTable % DV % Offsets(iP_bary)
  OS_S = EOSBaryonTable % DV % Offsets(iS_bary)
  OS_E = EOSBaryonTable % DV % Offsets(iE_bary)

  ALLOCATE( Ps_bary(SizeDs, SizeTs, SizeYps) )
  ALLOCATE( Ss_bary(SizeDs, SizeTs, SizeYps) )
  ALLOCATE( Es_bary(SizeDs, SizeTs, SizeYps) )

  Ps_bary = EOSBaryonTable % DV % Variables(iP_bary ) % Values
  Ss_bary = EOSBaryonTable % DV % Variables(iS_bary ) % Values
  Es_bary = EOSBaryonTable % DV % Variables(iE_bary ) % Values

  ALLOCATE( Ds_full(SizeDs) )
  ALLOCATE( Ts_full(SizeTs) )
  ALLOCATE( Yps_full(SizeYps) )
  
  Ds_full = Ds_bary
  Ts_full = Ts_bary
  Yps_full = Yps_bary
  
  ALLOCATE( Ps_full(SizeDs, SizeTs, SizeYps) )
  ALLOCATE( Ss_full(SizeDs, SizeTs, SizeYps) )
  ALLOCATE( Es_full(SizeDs, SizeTs, SizeYps) )
           
  ! Build full EOS for comparison purposes
  DO iD=1,SizeDs
    DO iT=1,SizeTs
      DO iYp=1,SizeYps

        Ymu_temp = Yps_bary(iYp) / Yp_over_Ymu
        Ymu_temp = 0.0d0

        ! Now add electron component
        ! Initialize temperature, DensitY, Ye
        ElectronPhotonState % t = Ts_bary(iT)
        ElectronPhotonState % rho = Ds_bary(iD)
        ElectronPhotonState % ye = Yps_bary(iYp) - Ymu_temp
        
        ! calculate electron quantities
        CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

        Es_helm = ElectronPhotonState % e
        Ps_helm = ElectronPhotonState % p
        Ss_helm = ElectronPhotonState % s

        ! calculate muon quantities 
        MuonState % t = Ts_bary(iT)
        MuonState % rhoym = Ds_bary(iD) * Ymu_temp
        
        CALL FullMuonEOS(MuonTable, MuonState)

        Es_muon = MuonState % e
        Ps_muon = MuonState % p
        Ss_muon = MuonState % s
                
        Es_full(iD,iT,iYp) = LOG10( 10.0d0**( Es_bary(iD,iT,iYp) ) + Es_helm + Es_muon)
        Ps_full(iD,iT,iYp) = LOG10( Ps_bary(iD,iT,iYp) + Ps_helm + Ps_muon)
        Ss_full(iD,iT,iYp) = LOG10( 10.0d0**( Ss_bary(iD,iT,iYp) ) + Ss_helm + Ss_muon)

        IF (Es_full(iD,iT,iYp) .NE. Es_full(iD,iT,iYp)) THEN
            WRITE(*,*) 'E', iD,iT,iYp, Es_bary(iD,iT,iYp), Es_helm, Es_muon
          STOP
        ENDIF
        IF (Ps_full(iD,iT,iYp) .NE. Ps_full(iD,iT,iYp)) THEN
          WRITE(*,*) 'P', iD,iT,iYp, Ps_bary(iD,iT,iYp), Ps_helm, Ps_muon
          STOP
        ENDIF

        IF (10.0d0**( Ps_full(iD,iT,iYp) ) - OS_P .NE. 10.0d0**( Ps_full(iD,iT,iYp)) - OS_P) THEN
          WRITE(*,*) '10**P', iD,iT,iYp, Ps_bary(iD,iT,iYp), Ps_helm, Ps_muon
          STOP
        ENDIF

        IF (Ss_full(iD,iT,iYp) .NE. Ss_full(iD,iT,iYp)) THEN
          WRITE(*,*) 'S', iD,iT,iYp, Ss_bary(iD,iT,iYp), Ss_helm, Ss_muon
          STOP
        ENDIF
                
      ENDDO
    ENDDO
  ENDDO

  CALL InitializeEOSComponentsInversion &
         ( Ds_bary, Ts_bary, Yps_bary, &
           10.0d0**( Es_full ) - OS_E, &
           10.0d0**( Ps_full ) - OS_P, &
           10.0d0**( Ss_full ) - OS_S, &
           BaryonPlusHelmTableName, &
           BaryonPlusHelmTableName, &
           Verbose_Option = .TRUE. )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET ENTER DATA &
  !$OMP MAP( to: Ds_bary, Ts_bary, Yps_bary, Es_bary, Ps_bary, Ss_bary, OS_E, OS_P, OS_S ) &
  !$OMP MAP( alloc: D, T, Yp, P, E, S, T_P, T_E, T_S, Error_P, Error_E, Error_S )
#elif defined (WEAKLIB_OACC)
  !$ACC ENTER DATA &
  !$ACC COPYIN( Ds_bary, Ts_bary, Yps_bary, Es_bary, Ps_bary, Ss_bary, OS_E, OS_P, OS_S ) &
  !$ACC CREATE( D, T, Yp, P, E, S, T_P, T_E, T_S, Error_P, Error_E, Error_S )
#endif

  WRITE(*,*)
  WRITE(*,'(A4,A10,I10.10)') '', 'nPoints = ', nPoints

  ! --- Random Sampling of EOSBaryonTable Table ---

  n_rndm = nPoints

  ! --- Initialize DensitY Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_D )

  associate &
    ( minD => EOSBaryonTable % TS % MinValues(iD_bary), &
      maxD => EOSBaryonTable % TS % MaxValues(iD_bary) )
  
  D(:) = 10.0d0**( LOG10(minD) + ( LOG10(maxD) - LOG10(minD) ) * rndm_D )

  PRINT*, "Min/Max D =  ", MINVAL( D ), MAXVAL( D )
  PRINT*, "             ", minD, maxD

  end associate

  ! --- Initialize Temperature Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_T )

  associate &
    ( minT => EOSBaryonTable % TS % MinValues(iT_bary), &
      maxT => EOSBaryonTable % TS % MaxValues(iT_bary) )

  T(:) = 10.0d0**( LOG10(minT) + ( LOG10(maxT) - LOG10(minT) ) * rndm_T )

  PRINT*, "Min/Max T  = ", MINVAL( T ), MAXVAL( T )
  PRINT*, "             ", minT, maxT

  end associate
  
  ! --- Initialize Electron Fraction Points ---

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_Yp )

  associate &
    ( minYp => EOSBaryonTable % TS % MinValues(iYp_bary), &
      maxYp => EOSBaryonTable % TS % MaxValues(iYp_bary) )
  
  Yp(:) = minYp + ( maxYp - minYp ) * rndm_Yp

  PRINT*, "Min/Max Yp = ", MINVAL( Yp ), MAXVAL( Yp )
  PRINT*, "             ", minYp, maxYp

  end associate

  ! --- Compute Internal Energy, Pressure, and Entropy ---

  DO iP = 1, nPoints
    
    Ym(iP) = Yp(iP) / Yp_over_Ymu
    Ym(iP) = 0.0d0
    Ye(iP) = Yp(iP) - Ym(iP)

    ! --- Interpolation of combined table ---
    CALL LogInterpolateSingleVariable_3D_Custom_Point &
           ( D(iP), T(iP), Yp(iP), Ds_full, Ts_full, Yps_full, OS_E, Es_full, E(iP) )

    CALL LogInterpolateSingleVariable_3D_Custom_Point &
           ( D(iP), T(iP), Yp(iP), Ds_full, Ts_full, Yps_full, OS_P, Ps_full, P(iP) )

    CALL LogInterpolateSingleVariable_3D_Custom_Point &
           ( D(iP), T(iP), Yp(iP), Ds_full, Ts_full, Yps_full, OS_S, Ss_full, S(iP) )
           
    ! --- Interpolation of separate table ---
    CALL LogInterpolateSingleVariable_3D_Custom_Point &
           ( D(iP), T(iP), Yp(iP), Ds_bary, Ts_bary, Yps_bary, OS_E, Es_bary, E_bary )

    iD  = Index1D_Log(  D(iP),  Ds_bary )
    iT  = Index1D_Log(  T(iP),  Ts_bary )
    iYp = Index1D_Lin( Yp(iP), Yps_bary )

    LocalOffset = MINVAL( Ps_bary(iD:iD+1,iT:iT+1,iYp:iYp+1) )
    IF (LocalOffset .lt. 0.0_dp) THEN
        LocalOffset = -1.1d0*LocalOffset
    ELSE
        LocalOffset = 0.0_dp
    ENDIF
    CALL LogInterpolateSingleVariable_3D_Custom_Point &
          ( D(iP)           , T(iP)           , Yp(iP)             , &
            Ds_bary(iD:iD+1), Ts_bary(iT:iT+1), Yps_bary(iYp:iYp+1), &
            LocalOffset, LOG10(Ps_bary(iD:iD+1,iT:iT+1,iYp:iYp+1) + LocalOffset), P_bary )

    CALL LogInterpolateSingleVariable_3D_Custom_Point &
           ( D(iP), T(iP), Yp(iP), Ds_bary, Ts_bary, Yps_bary, OS_S, Ss_bary, S_bary )
           
    ! Now add electron component
    ! Initialize temperature, DensitY, Ye
    ElectronPhotonState % t = T(iP)
    ElectronPhotonState % rho = D(iP)
    ElectronPhotonState % ye = Ye(iP)
    
    ! calculate electron quantities
    CALL ElectronPhotonEOS(HelmholtzTable, ElectronPhotonState)

    Es_helm = ElectronPhotonState % e
    Ps_helm = ElectronPhotonState % p
    Ss_helm = ElectronPhotonState % s

    ! calculate muon quantities 
    MuonState % t = T(iP)
    MuonState % rhoym = D(iP) * Ym(iP)
    
    CALL FullMuonEOS(MuonTable, MuonState)

    Es_muon = MuonState % e
    Ps_muon = MuonState % p
    Ss_muon = MuonState % s
           
    E_int_error(iP) = ABS( E_bary + Es_helm + Es_muon - E(iP))/ABS(E(iP))
    P_int_error(iP) = ABS( P_bary + Ps_helm + Ps_muon - P(iP))/ABS(P(iP))
    S_int_error(iP) = ABS( S_bary + Ss_helm + Ss_muon - S(iP))/ABS(S(iP))

#ifdef INTERPOLATION_SPLIT_TABLE_COMBINED
#else
    E(iP) = E_bary + Es_helm + Es_muon
    P(iP) = P_bary + Ps_helm + Ps_muon
    S(iP) = S_bary + Ss_helm + Ss_muon
#endif

  END DO
  
  WRITE(*,*)
  WRITE(*,*) "Min/Max E = ", MINVAL( E ), MAXVAL( E )
  WRITE(*,*) "Min/Max P = ", MINVAL( P ), MAXVAL( P )
  WRITE(*,*) "Min/Max S = ", MINVAL( S ), MAXVAL( S )

  WRITE(*,*)
  WRITE(*,*) 'Interpolation errors for PES'
  WRITE(*,*) MAXVAL(P_int_error), MINVAL(P_int_error)
  WRITE(*,*) MAXVAL(E_int_error), MINVAL(E_int_error)
  WRITE(*,*) MAXVAL(S_int_error), MINVAL(S_int_error)

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET UPDATE TO( D, T, Yp, P, E, S )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC UPDATE DEVICE( D, T, Yp, P, E, S )
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
    ( minT => 1.0001 * EOSBaryonTable % TS % MinValues(iT_bary), &
      maxT => 0.9999 * EOSBaryonTable % TS % MaxValues(iT_bary) )

  T_E = MIN( MAX( T + dT, minT ), maxT )

  end associate

  Error_E = 0

  CALL CPU_TIME( tBegin )

! #if defined(WEAKLIB_OMP_OL)
  ! !$OMP TARGET UPDATE TO( T_E, Error_E )
! #elif defined (WEAKLIB_OACC)
  ! !$ACC UPDATE DEVICE( T_E, Error_E )
! #endif

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  !$OMP PRIVATE( T_Guess )
#elif defined (WEAKLIB_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR &
  !$ACC PRIVATE( T_Guess ) &
  !$ACC PRESENT( D, E, Yp, Ds_bary, Ts_bary, Yps_bary, Es_bary, OS_E, T_E, Error_E )
#elif defined (WEAKLIB_OMP)
  !$OMP PARALLEL DO &
  !$OMP PRIVATE( T_Guess )
#endif
  DO iP = 1, nPoints
    T_Guess = T_E(iP)
    CALL ComputeTemperatureWith_DEYpYl_Single_Guess &
           ( D(iP), E(iP), Ye(iP), Ym(iP), &
           Ds_bary, Ts_bary, Yps_bary, Es_bary, OS_E, &
           T_E(iP), T_Guess, Error_E(iP) )
  END DO
  
  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  PRINT*
  PRINT*, "ComputeTemperatureWith_DEYp (Good Guess):"
  PRINT*, "CPU_TIME = ", tCPU

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE FROM( T_E, Error_E )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE HOST( T_E, Error_E )
#endif
  DO iP = 1, nPoints
    IF( Error_E(iP) .NE. 0 ) THEN 
      WRITE(*,*) iP, D(iP), E(iP), T(iP), Ye(iP), Ym(iP)
      CALL DescribeEOSComponentsInversionError( Error_E(iP) )
    ENDIF
  END DO

  Error = ABS( T - T_E ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_E(iMaxError), Yp(iMaxError), Ds_full, Ts_full, Yps_full, OS_E, Es_full, E_T )
  PRINT*, "CPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "CPU  Error(E) = ", ABS( E(iMaxError) - E_T ) / E(iMaxError)
  PRINT*, "CPU MINVAL(T) = ", MINVAL( T_E )

  ! -------------------------------------------------------------------

  associate &
    ( minT => 1.0001 * EOSBaryonTable % TS % MinValues(iT_bary), &
      maxT => 0.9999 * EOSBaryonTable % TS % MaxValues(iT_bary) )

  Amp = 0.1_dp * ( maxT - minT ) ! --- Perturbation Amplitude

  PRINT*
  PRINT*, "Perturbation Amplitude: ", Amp

  dT = Amp * 2.0_dp * ( rndm_T - 0.5_dp ) * T

  T_E = MIN( MAX( T + dT, minT ), maxT )

  end associate

  Error_E = 0

  CALL CPU_TIME( tBegin )


#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  !$OMP PRIVATE( T_Guess )
#elif defined (WEAKLIB_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR &
  !$ACC PRIVATE( T_Guess ) &
  !$ACC PRESENT( D, E, Yp, Ds_bary, Ts_bary, Yps_bary, Es_bary, OS_E, T_E, Error_E )
#elif defined (WEAKLIB_OMP)
  !$OMP PARALLEL DO &
  !$OMP PRIVATE( T_Guess )
#endif
  DO iP = 1, nPoints
    T_Guess = T_E(iP)
    CALL ComputeTemperatureWith_DEYpYl_Single_Guess &
           ( D(iP), E(iP), Ye(iP), Ym(iP), &
             Ds_bary, Ts_bary, Yps_bary, Es_bary, OS_E, T_E(iP), T_Guess, &
             Error_E(iP) )
  END DO

  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  PRINT*
  PRINT*, "ComputeTemperatureWith_DEYp (Bad Guess):"
  PRINT*, "CPU_TIME = ", tCPU

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE FROM( T_E, Error_E )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE HOST( T_E, Error_E )
#endif
  DO iP = 1, nPoints
    IF( Error_E(iP) .NE. 0 ) CALL DescribeEOSComponentsInversionError( Error_E(iP) )
  END DO
  Error = ABS( T - T_E ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_E(iMaxError), Yp(iMaxError), Ds_full, Ts_full, Yps_full, OS_E, Es_full, E_T )
  PRINT*, "CPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "CPU  Error(E) = ", ABS( E(iMaxError) - E_T ) / E(iMaxError)
  PRINT*, "CPU MINVAL(T) = ", MINVAL( T_E )

  ! -------------------------------------------------------------------

  T_E = 0.0_dp
  Error_E = 0

  CALL CPU_TIME( tBegin )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  !$OMP PRIVATE( T_Guess )
#elif defined (WEAKLIB_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR &
  !$ACC PRIVATE( T_Guess ) &
  !$ACC PRESENT( D, E, Yp, Ds_bary, Ts_bary, Yps_bary, Es_bary, OS_E, T_E, Error_E )
#elif defined (WEAKLIB_OMP)
  !$OMP PARALLEL DO &
  !$OMP PRIVATE( T_Guess )
#endif
  DO iP = 1, nPoints
    T_Guess = T_E(iP)
    CALL ComputeTemperatureWith_DEYpYl_Single_NoGuess &
           ( D(iP), E(iP), Ye(iP), Ym(iP), &
             Ds_bary, Ts_bary, Yps_bary, Es_bary, OS_E, T_E(iP), &
             Error_E(iP) )
  END DO


  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  PRINT*
  PRINT*, "ComputeTemperatureWith_DEYp (No Guess):"
  PRINT*, "CPU_TIME = ", tCPU

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE FROM( T_E, Error_E )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE HOST( T_E, Error_E )
#endif
  DO iP = 1, nPoints
    IF( Error_E(iP) .NE. 0 ) CALL DescribeEOSComponentsInversionError( Error_E(iP) )
  END DO
  Error = ABS( T - T_E ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_E(iMaxError), Yp(iMaxError), Ds_full, Ts_full, Yps_full, OS_E, Es_full, E_T )
  PRINT*, "CPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "CPU  Error(E) = ", ABS( E(iMaxError) - E_T ) / E(iMaxError)
  PRINT*, "CPU MINVAL(T) = ", MINVAL( T_E )
  
  ! -------------------------------------------------------------------
  ! --- Recover Temperature from Entropy ----------------------
  ! -------------------------------------------------------------------

  Amp = 0.001_dp ! --- Perturbation Amplitude
  dT = Amp * 2.0_dp * ( rndm_T - 0.5_dp ) * T

  PRINT*, "Perturbation Amplitude: ", Amp

  associate &
    ( minT => 1.0001 * EOSBaryonTable % TS % MinValues(iT_bary), &
      maxT => 0.9999 * EOSBaryonTable % TS % MaxValues(iT_bary) )

  T_S = MIN( MAX( T + dT, minT ), maxT )

  end associate

  Error_S = 0

  CALL CPU_TIME( tBegin )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  !$OMP PRIVATE( T_Guess )
#elif defined (WEAKLIB_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR &
  !$ACC PRIVATE( T_Guess ) &
  !$ACC PRESENT( D, P, Yp, Ds_bary, Ts_bary, Yps_bary, Ss_bary, OS_S, T_S, Error_S )
#elif defined(WEAKLIB_OMP_OL)
  !$OMP PARALLEL DO &
  !$OMP PRIVATE( T_Guess )
#endif
  DO iP = 1, nPoints
    T_Guess = T_S(iP)
    CALL ComputeTemperatureWith_DSYpYl_Single_Guess &
           ( D(iP), S(iP), Ye(iP), Ym(iP), &
             Ds_bary, Ts_bary, Yps_bary, Ss_bary, OS_S, T_S(iP), T_Guess, &
             Error_S(iP) )
  END DO
  WRITE(*,*) 1, MAXVAL(Error_S), MAXVAL(T_S)

  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  PRINT*
  PRINT*, "ComputeTemperatureWith_DSYp (Good Guess):"
  PRINT*, "CPU_TIME = ", tCPU

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE FROM( T_S, Error_S )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE HOST( T_S, Error_S )
#endif
  DO iP = 1, nPoints
    IF( Error_S(iP) .NE. 0 ) CALL DescribeEOSComponentsInversionError( Error_S(iP) )
  END DO
  
  Error = ABS( T - T_S ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_S(iMaxError), Yp(iMaxError), Ds_full, Ts_full, Yps_full, OS_S, Ss_full, S_T )
  PRINT*, "CPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "CPU  Error(S) = ", ABS( S(iMaxError) - S_T ) / S(iMaxError)
  PRINT*, "CPU MINVAL(T) = ", MINVAL( T_S )

  ! -------------------------------------------------------------------

  associate &
    ( minT => 1.0001 * EOSBaryonTable % TS % MinValues(iT_bary), &
      maxT => 0.9999 * EOSBaryonTable % TS % MaxValues(iT_bary) )

  Amp = 0.1_dp * ( maxT - minT )
  dT = Amp * 2.0_dp * ( rndm_T - 0.5_dp ) * T

  PRINT*, "Perturbation Amplitude: ", Amp

  T_S = MIN( MAX( T + dT, minT ), maxT )

  end associate

  Error_S = 0

  CALL CPU_TIME( tBegin )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  !$OMP PRIVATE( T_Guess )
#elif defined (WEAKLIB_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR &
  !$ACC PRIVATE( T_Guess ) &
  !$ACC PRESENT( D, S, Yp, Ds_bary, Ts_bary, Yps_bary, Ss_bary, OS_S, T_S, Error_S )
#elif defined (WEAKLIB_OMP)
  !$OMP PARALLEL DO &
  !$OMP PRIVATE( T_Guess )
#endif
  DO iP = 1, nPoints
    T_Guess = T_S(iP)
    CALL ComputeTemperatureWith_DSYpYl_Single_Guess &
           ( D(iP), S(iP), Ye(iP), Ym(iP), &
             Ds_bary, Ts_bary, Yps_bary, Ss_bary, OS_S, T_S(iP), T_Guess, &
             Error_S(iP) )
  END DO
  
  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  PRINT*
  PRINT*, "ComputeTemperatureWith_DSYp (Bad Guess):"
  PRINT*, "CPU_TIME = ", tCPU

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE FROM( T_S, Error_S )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE HOST( T_S, Error_S )
#endif
  DO iP = 1, nPoints
    IF( Error_S(iP) .NE. 0 ) CALL DescribeEOSComponentsInversionError( Error_S(iP) )
  END DO
  
  Error = ABS( T - T_S ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_S(iMaxError), Yp(iMaxError), Ds_full, Ts_full, Yps_full, OS_S, Ss_full, S_T )
  PRINT*, "CPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "CPU  Error(S) = ", ABS( S(iMaxError) - S_T ) / S(iMaxError)
  PRINT*, "CPU MINVAL(T) = ", MINVAL( T_S )

  ! -------------------------------------------------------------------

  T_S = 0.0_dp
  Error_S = 0

  CALL CPU_TIME( tBegin )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  !$OMP PRIVATE( T_Guess )
#elif defined (WEAKLIB_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR &
  !$ACC PRIVATE( T_Guess ) &
  !$ACC PRESENT( D, S, Yp, Ds_bary, Ts_bary, Yps_bary, Ss_bary, OS_S, T_S, Error_S )
#elif defined (WEAKLIB_OMP)
  !$OMP PARALLEL DO &
  !$OMP PRIVATE( T_Guess )
#endif
  DO iP = 1, nPoints
    T_Guess = T_S(iP)
    CALL ComputeTemperatureWith_DSYpYl_Single_NoGuess &
           ( D(iP), S(iP), Ye(iP), Ym(iP), &
             Ds_bary, Ts_bary, Yps_bary, Ss_bary, OS_S, T_S(iP), &
             Error_S(iP) )
  END DO

  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  PRINT*
  PRINT*, "ComputeTemperatureWith_DSYp (No Guess):"
  PRINT*, "CPU_TIME = ", tCPU

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE FROM( T_S, Error_S )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE HOST( T_S, Error_S )
#endif
  DO iP = 1, nPoints
    IF( Error_S(iP) .NE. 0 ) CALL DescribeEOSComponentsInversionError( Error_S(iP) )
  END DO
  
  Error = ABS( T - T_S ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_S(iMaxError), Yp(iMaxError), Ds_full, Ts_full, Yps_full, OS_S, Ss_full, S_T )
  PRINT*, "CPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "CPU  Error(S) = ", ABS( S(iMaxError) - S_T ) / S(iMaxError)
  PRINT*, "CPU MINVAL(T) = ", MINVAL( T_S )

  ! -------------------------------------------------------------------
  ! --- Recover Temperature from Pressure ----------------------
  ! -------------------------------------------------------------------

  Amp = 0.001_dp ! --- Perturbation Amplitude

  PRINT*
  PRINT*, "Perturbation Amplitude: ", Amp

  CALL RANDOM_SEED( SIZE = n_rndm )
  CALL RANDOM_NUMBER( rndm_T )
  
  dT = Amp * 2.0_dp * ( rndm_T - 0.5_dp ) * T

  associate &
    ( minT => 1.0001 * EOSBaryonTable % TS % MinValues(iT_bary), &
      maxT => 0.9999 * EOSBaryonTable % TS % MaxValues(iT_bary) )

  T_P = MIN( MAX( T + dT, minT ), maxT )

  end associate

  Error_P = 0

  CALL CPU_TIME( tBegin )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  !$OMP PRIVATE( T_Guess )
#elif defined (WEAKLIB_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR &
  !$ACC PRIVATE( T_Guess ) &
  !$ACC PRESENT( D, P, Yp, Ds_bary, Ts_bary, Yps_bary, Ps_bary, OS_P, T_P, Error_P )
#elif defined(WEAKLIB_OMP_OL)
  !$OMP PARALLEL DO &
  !$OMP PRIVATE( T_Guess )
#endif
  DO iP = 1, nPoints
    T_Guess = T_P(iP)
    CALL ComputeTemperatureWith_DPYpYl_Single_Guess &
           ( D(iP), P(iP), Ye(iP), Ym(iP), &
             Ds_bary, Ts_bary, Yps_bary, Ps_bary, OS_P, T_P(iP), T_Guess, &
             Error_P(iP) )
  END DO

  CALL CPU_TIME( tEnd )
  
  tCPU = tEnd - tBegin

  PRINT*
  PRINT*, "ComputeTemperatureWith_DPYp (Good Guess):"
  PRINT*, "CPU_TIME = ", tCPU

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE FROM( T_P, Error_P )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE HOST( T_P, Error_P )
#endif
  DO iP = 1, nPoints
    IF( Error_P(iP) .NE. 0 ) CALL DescribeEOSComponentsInversionError( Error_P(iP) )
  END DO
  
  Error = ABS( T - T_P ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_P(iMaxError), Yp(iMaxError), Ds_full, Ts_full, Yps_full, OS_P, Ps_full, P_T )
  PRINT*, "CPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "CPU  Error(P) = ", ABS( P(iMaxError) - P_T ) / P(iMaxError)
  PRINT*, "CPU MINVAL(T) = ", MINVAL( T_P )

  ! -------------------------------------------------------------------

  associate &
    ( minT => 1.0001 * EOSBaryonTable % TS % MinValues(iT_bary), &
      maxT => 0.9999 * EOSBaryonTable % TS % MaxValues(iT_bary) )

  Amp = 0.1_dp * ( maxT - minT ) ! --- Perturbation Amplitude

  PRINT*
  PRINT*, "Perturbation Amplitude: ", Amp

  dT = Amp * 2.0_dp * ( rndm_T - 0.5_dp ) * T

  T_P = MIN( MAX( T + dT, minT ), maxT )

  end associate

  Error_P = 0

  CALL CPU_TIME( tBegin )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  !$OMP PRIVATE( T_Guess )
#elif defined (WEAKLIB_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR &
  !$ACC PRIVATE( T_Guess ) &
  !$ACC PRESENT( D, P, Yp, Ds_bary, Ts_bary, Yps_bary, Ps_bary, OS_P, T_P, Error_P )
#elif defined(WEAKLIB_OMP_OL)
  !$OMP PARALLEL DO &
  !$OMP PRIVATE( T_Guess )
#endif
  DO iP = 1, nPoints
    T_Guess = T_P(iP)
    CALL ComputeTemperatureWith_DPYpYl_Single_Guess &
           ( D(iP), P(iP), Ye(iP), Ym(iP), &
             Ds_bary, Ts_bary, Yps_bary, Ps_bary, OS_P, T_P(iP), T_Guess, &
             Error_P(iP) )
  END DO

  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  PRINT*
  PRINT*, "ComputeTemperatureWith_DPYp (Bad Guess):"
  PRINT*, "CPU_TIME = ", tCPU

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE FROM( T_P, Error_P )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE HOST( T_P, Error_P )
#endif
  DO iP = 1, nPoints
    IF( Error_P(iP) .NE. 0 ) CALL DescribeEOSComponentsInversionError( Error_P(iP) )
  END DO

  Error = ABS( T - T_P ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_P(iMaxError), Yp(iMaxError), Ds_full, Ts_full, Yps_full, OS_P, Ps_full, P_T )
  PRINT*, "CPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "CPU  Error(P) = ", ABS( P(iMaxError) - P_T ) / P(iMaxError)
  PRINT*, "CPU MINVAL(T) = ", MINVAL( T_P )

  ! -------------------------------------------------------------------

  T_P = 0.0_dp
  Error_P = 0

  CALL CPU_TIME( tBegin )

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD &
  !$OMP PRIVATE( T_Guess )
#elif defined (WEAKLIB_OACC)
  !$ACC PARALLEL LOOP GANG VECTOR &
  !$ACC PRIVATE( T_Guess ) &
  !$ACC PRESENT( D, P, Yp, Ds_bary, Ts_bary, Yps_bary, Ps_bary, OS_P, T_P, Error_P )
#elif defined (WEAKLIB_OMP)
  !$OMP PARALLEL DO &
  !$OMP PRIVATE( T_Guess )
#endif
  DO iP = 1, nPoints
    T_Guess = T_P(iP)
    CALL ComputeTemperatureWith_DPYpYl_Single_NoGuess &
           ( D(iP), P(iP), Ye(iP), Ym(iP), &
             Ds_bary, Ts_bary, Yps_bary, Ps_bary, OS_P, T_P(iP), &
             Error_P(iP) )
  END DO

  CALL CPU_TIME( tEnd )
  tCPU = tEnd - tBegin

  PRINT*
  PRINT*, "ComputeTemperatureWith_DPYp (No Guess):"
  PRINT*, "CPU_TIME = ", tCPU

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET UPDATE FROM( T_P, Error_P )
#elif defined (WEAKLIB_OACC)
  !$ACC UPDATE HOST( T_P, Error_P )
#endif
  DO iP = 1, nPoints
    IF( Error_P(iP) .NE. 0 ) CALL DescribeEOSComponentsInversionError( Error_P(iP) )
  END DO

  Error = ABS( T - T_P ) / T
  iMaxError = MAXLOC( Error, DIM=1 )
  CALL LogInterpolateSingleVariable_3D_Custom_Point &
         ( D(iMaxError), T_P(iMaxError), Yp(iMaxError), Ds_full, Ts_full, Yps_full, OS_P, Ps_full, P_T )
  PRINT*, "CPU  Error(T) = ", MAXVAL( Error ), MINVAL( Error ), SUM( Error ) / DBLE( nPoints )
  PRINT*, "CPU  Error(P) = ", ABS( P(iMaxError) - P_T ) / P(iMaxError)
  PRINT*, "CPU MINVAL(T) = ", MINVAL( T_P )

  ! -------------------------------------------------------------------

#if defined(WEAKLIB_OMP_OL)
  !$OMP TARGET EXIT DATA &
  !$OMP MAP( release: D, E, P, S, Yp, Ds_bary, Ts_bary, Yps_bary, Es_bary, Ps_bary, Ss_bary, OS_E, OS_P, OS_S, T_S, Error_E )
#elif defined (WEAKLIB_OACC)
  !$ACC EXIT DATA &
  !$ACC DELETE( D, E, P, S, Yp, Ds_bary, Ts_bary, Yps_bary, Es_bary, Ps_bary, Ss_bary, OS_E, OS_P, OS_S, T_S, Error_E )
#endif

  DEALLOCATE( Ds_bary, Ts_bary, Yps_bary, Ps_bary, Ss_bary, Es_bary )
  DEALLOCATE( Ds_full, Ts_full, Yps_full, Ps_full, Ss_full, Es_full )

#if defined(WEAKLIB_OMP_OL)
#elif defined(WEAKLIB_OACC)
  !$ACC SHUTDOWN
#endif

  CALL MPI_FINALIZE( ierr )

END PROGRAM wlComposeNewInversionTest
